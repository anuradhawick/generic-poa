use crate::alignment::SeqGraphAlignment;
use petgraph::{
    algo::toposort,
    graph::{DiGraph, NodeIndex},
    visit::EdgeRef,
    Direction,
};
use std::cmp::max;

#[derive(Debug)]
pub struct NodeData {
    pub item: String,
    pub aligned_to: Vec<NodeIndex>,
}

#[derive(Debug)]
pub struct EdgeData {
    pub labels: Vec<String>,
}

pub struct POAGraph {
    pub graph: DiGraph<NodeData, EdgeData>,
    pub sequeces: Vec<Vec<String>>,
    pub labels: Vec<String>,
    pub start_indices: Vec<NodeIndex>,
    pub width: usize,
}

impl POAGraph {
    /// Initialise the POA graph with the first sequence
    pub fn new(label: String, seq: Vec<String>) -> Self {
        let mut width = 0usize;
        let mut graph = DiGraph::new();
        let nodes: Vec<NodeIndex> = seq
            .iter()
            .map(|item| {
                width = max(item.len(), width);
                graph.add_node(NodeData {
                    item: item.clone(),
                    aligned_to: vec![],
                })
            })
            .collect();
        for (position, &index) in nodes.iter().enumerate() {
            if position < nodes.len() - 1 {
                graph.add_edge(
                    index,
                    nodes[position + 1],
                    EdgeData {
                        labels: vec![label.clone()],
                    },
                );
            }
        }
        Self {
            graph,
            sequeces: vec![seq],
            labels: vec![label],
            start_indices: vec![nodes[0]],
            width,
        }
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

        // println!("valid_seq_positions={valid_seq_positions:?}");
        // println!("seq_start_pos={seq_start_pos:?}");
        // println!("seq_end_pos={seq_end_pos:?}");

        if seq_start_pos > 0 {
            (first_node_index, head_node_index) =
                self.add_seq_segment(aln.label.clone(), &seq[0..seq_start_pos as usize]);
        }
        if seq_end_pos < seq.len() as i32 {
            (tail_node_index, _) =
                self.add_seq_segment(aln.label.clone(), &seq[seq_end_pos as usize + 1..]);
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
            // a gap
            if seq_pos.is_none() {
                continue;
            }

            let seq_item = &seq[seq_pos.unwrap() as usize];
            let node_index;

            // this is an aligned position
            if let Some(match_node_index) = match_node_index {
                // matching alignment
                if self.graph[match_node_index].item == *seq_item {
                    node_index = match_node_index;
                }
                // non matching alignment
                else {
                    let other_aligned = &self.graph[match_node_index].aligned_to;
                    let mut found_node = None;
                    for other_node in other_aligned {
                        if self.graph[*other_node].item == *seq_item {
                            found_node = Some(other_node);
                        }
                    }
                    // the mismatch is already accounted
                    if let Some(found_node) = found_node {
                        node_index = *found_node;
                    }
                    // a new mismatching base
                    else {
                        let other_node_indices =
                            [vec![match_node_index].as_slice(), other_aligned.as_slice()].concat();
                        for other_node_index in other_node_indices.iter() {
                            self.graph[*other_node_index]
                                .aligned_to
                                .push(*other_node_index);
                        }
                        node_index = self.graph.add_node(NodeData {
                            item: seq_item.clone(),
                            aligned_to: other_node_indices,
                        });
                    }
                }
            }
            // not aligned, this is a new base insertion
            else {
                node_index = self.graph.add_node(NodeData {
                    item: seq_item.clone(),
                    aligned_to: vec![],
                })
            }

            // if a new start is there
            if let Some(head_node_index) = head_node_index {
                self.add_or_update_edge(head_node_index, node_index, aln.label.clone())
            }

            // update head
            head_node_index = Some(node_index);

            // update first node
            if first_node_index.is_none() {
                first_node_index = head_node_index;
            }
        }

        // add the edges
        if let (Some(head_node_index), Some(tail_node_index)) = (head_node_index, tail_node_index) {
            self.add_or_update_edge(head_node_index, tail_node_index, aln.label.clone());
        }

        // record the summaries
        self.sequeces.push(seq);
        self.labels.push(aln.label);
        self.start_indices.push(first_node_index.unwrap());
    }

    pub fn get_string(&self) -> String {
        let indices = toposort(&self.graph, None).unwrap();
        let mut output = String::new();
        for index in indices {
            output.push_str(&format!("{:?}:{}\n", index.index(), self.graph[index].item));

            for edge in self.graph.edges_directed(index, Direction::Outgoing) {
                output.push_str(&format!(
                    "\t{:?} -> {:?} {:?}\n",
                    edge.source().index(),
                    edge.target().index(),
                    edge.weight().labels
                ));
            }
        }

        output
    }

    fn add_or_update_edge(&mut self, a: NodeIndex, b: NodeIndex, label: String) {
        if let Some(edge) = self.graph.find_edge(a, b) {
            self.graph[edge].labels.push(label);
        } else {
            self.graph.add_edge(
                a,
                b,
                EdgeData {
                    labels: vec![label],
                },
            );
        }
    }

    fn add_seq_segment(
        &mut self,
        label: String,
        seq: &[String],
    ) -> (Option<NodeIndex>, Option<NodeIndex>) {
        // println!("Adding segment label={label} seq={seq:?}");
        let nodes: Vec<NodeIndex> = seq
            .iter()
            .map(|item| {
                self.width = max(item.len(), self.width);
                self.graph.add_node(NodeData {
                    item: item.clone(),
                    aligned_to: vec![],
                })
            })
            .collect();

        for (position, &index) in nodes.iter().enumerate() {
            if position < nodes.len() - 1 {
                self.add_or_update_edge(index, nodes[position + 1], label.clone());
            }
        }
        (nodes.first().copied(), nodes.last().copied())
    }
}

#[cfg(test)]
mod graph_tests {
    use super::*;

    #[test]
    fn new_test() {
        let seq = vec!["ABC".to_string(), "BBC".to_string(), "DDD".to_string()];
        let graph = POAGraph::new("seq_1".to_string(), seq);
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
        let mut graph = POAGraph::new("seq_1".to_string(), seq1);
        let sg_aln = SeqGraphAlignment::align_seq_to_graph("seq_2".to_string(), seq2, &graph.graph);
        graph.add_alignment(sg_aln);
        println!("{}", graph.get_string());
    }
}
