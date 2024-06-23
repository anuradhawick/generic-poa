use crate::alignment::SeqGraphAlignment;
use ndarray::Array2;
use petgraph::{
    algo::toposort,
    data::Build,
    graph::{DiGraph, NodeIndex},
    visit::EdgeRef,
    Direction,
};
use std::cmp::max;
use std::collections::HashMap;

#[derive(Debug)]
pub struct NodeData {
    pub item: String,
    pub position: usize,
    pub aligned_to: Vec<NodeIndex>,
}

#[derive(Debug)]
pub struct EdgeData {}

pub struct POAGraph {
    pub graph: DiGraph<NodeData, EdgeData>,
}

impl POAGraph {
    /// Initialise the POA graph with the first sequence
    pub fn new(seq: Vec<String>) -> Self {
        let mut graph = DiGraph::new();
        let nodes: Vec<NodeIndex> = seq
            .into_iter()
            .enumerate()
            .map(|(position, item)| {
                graph.add_node(NodeData {
                    item,
                    position,
                    aligned_to: vec![],
                })
            })
            .collect();
        for (position, &index) in nodes.iter().enumerate() {
            if position < nodes.len() - 1 {
                graph.add_edge(index, nodes[position + 1], EdgeData {});
            }
        }
        Self { graph }
    }

    pub fn add_alignment(&mut self, aln: SeqGraphAlignment) {
        let seq = aln.seq;
        let seq_match_positions = aln.seq_match_positions;
        let graph_match_node_indices = aln.graph_match_node_indices;

        let mut first_node_index = None;
        let mut head_node_index = None;
        let mut tail_node_index = None;

        let valid_seq_positions: Vec<i32> = seq_match_positions
            .iter()
            .filter_map(|x| x.as_ref().map(|n| *n))
            .collect();
        let seq_start_pos = *valid_seq_positions.first().unwrap();
        let seq_end_pos = *valid_seq_positions.last().unwrap();

        if seq_start_pos > 0 {
            let first_and_head = self.add_seq_segment(&seq[0..seq_start_pos as usize]);
            first_node_index = Some(first_and_head.0);
            head_node_index = Some(first_and_head.1);
        }
        if seq_end_pos < seq.len() as i32 {
            let tail_segment = self.add_seq_segment(&seq[seq_end_pos as usize + 1..]);
            tail_node_index = Some(tail_segment.0);
        }

        // now we march along the aligned part. For each base, we find or create
        // a node in the graph:
        //   - if unmatched, the corresponding node is a new node
        //   - if matched:
        //       - if matched to a node with the same base, the node is that node
        //       - if matched to a node with a different base whch is in turn
        //         aligned to a node with the same base, that aligned node is
        //         the node
        //       - otherwise, we create a new node.
        // In all cases, we create edges (or add labels) threading through the
        // nodes.

        for (&seq_pos, &match_node_index) in
            seq_match_positions.iter().zip(&graph_match_node_indices)
        {
            if seq_pos.is_none() {
                continue;
            }

            let seq_item = &seq[seq_pos.unwrap() as usize];
            let node_index;
            if let Some(match_node_index) = match_node_index {
                if self.graph[match_node_index].item == *seq_item {
                    node_index = match_node_index;
                } else {
                    let other_aligned = &self.graph[match_node_index].aligned_to;
                    let mut found_node = None;
                    for other_node in other_aligned {
                        if self.graph[*other_node].item == *seq_item {
                            found_node = Some(other_node);
                        }
                    }
                    if let Some(found_node) = found_node {
                        node_index = *found_node;
                    } else {
                        let other_node_indices =
                            [vec![match_node_index].as_slice(), other_aligned.as_slice()].concat();
                        for other_node_index in other_node_indices.iter() {
                            self.graph[*other_node_index]
                                .aligned_to
                                .push(*other_node_index);
                        }
                        node_index = self.graph.add_node(NodeData {
                            item: seq_item.clone(),
                            position: 0,
                            aligned_to: other_node_indices,
                        });
                    }
                }
            } else {
                node_index = self.graph.add_node(NodeData {
                    item: seq_item.clone(),
                    position: 0,
                    aligned_to: vec![],
                })
            }
            if let Some(head_node_index) = head_node_index {
                self.graph
                    .add_edge(head_node_index, node_index, EdgeData {});
            }
            head_node_index = Some(node_index);

            if first_node_index.is_none() {
                first_node_index = head_node_index;
            }
        }

        if let (Some(head_node_index), Some(tail_node_index)) = (head_node_index, tail_node_index) {
            self.graph
                .add_edge(head_node_index, tail_node_index, EdgeData {});
        }
    }

    fn add_seq_segment(&mut self, seq: &[String]) -> (NodeIndex, NodeIndex) {
        let nodes: Vec<NodeIndex> = seq
            .iter()
            .enumerate()
            .map(|(position, item)| {
                self.graph.add_node(NodeData {
                    item: item.clone(),
                    position,
                    aligned_to: vec![],
                })
            })
            .collect();

        for (position, &index) in nodes.iter().enumerate() {
            if position < nodes.len() - 1 {
                self.graph.add_edge(index, nodes[position + 1], EdgeData {});
            }
        }

        let first = *nodes.first().unwrap();
        let last = *nodes.last().unwrap();
        (first, last)
    }
}

#[cfg(test)]
mod graph_tests {
    use super::*;

    #[test]
    fn new_test() {
        let seq = vec!["ABC".to_string(), "BBC".to_string(), "DDD".to_string()];
        let graph = POAGraph::new(seq);
        assert_eq!(graph.graph.node_indices().len(), 3);
        assert_eq!(graph.graph.edge_indices().len(), 2);
    }

    #[test]
    fn add_aln_test() {
        let seq1 = vec![
            "T".to_string(),
            "G".to_string(),
            "X".to_string(),
            "T".to_string(),
        ];
        let seq2 = vec![
            "A".to_string(),
            "T".to_string(),
            "G".to_string(),
            "X".to_string(),
            "T".to_string(),
        ];
        let mut graph = POAGraph::new(seq1);
        let sg_aln = SeqGraphAlignment::align_seq_to_graph(seq2, &graph.graph);
        graph.add_alignment(sg_aln);
    }
}
