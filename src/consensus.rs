use crate::graph::{EdgeData, NodeData};
use petgraph::{
    algo::toposort,
    graph::{DiGraph, NodeIndex},
    visit::EdgeRef,
    Direction,
};
use std::collections::HashMap;

pub struct Consensus {
    graph: DiGraph<NodeData, EdgeData>,
    start_indices: Vec<NodeIndex>,
    labels: Vec<String>,
}

impl Consensus {
    pub fn new(
        graph: DiGraph<NodeData, EdgeData>,
        start_indices: Vec<NodeIndex>,
        labels: Vec<String>,
    ) -> Self {
        Self {
            graph,
            start_indices,
            labels,
        }
    }

    pub fn compute(&self) -> Vec<(String, Vec<String>)> {
        // Step 1: assign node IDs to columns in the output
        //      column_index[node.ID] is the position in the toposorted node list
        //      of the node itself, or the earliest node it is aligned to.
        let indices = toposort(&self.graph, None).unwrap();
        let mut column_index: HashMap<NodeIndex, i32> = HashMap::new();
        let mut current_column = 0;

        for index in indices {
            let other_columns: Vec<i32> = self.graph[index]
                .aligned_to
                .iter()
                .filter_map(|other_index| column_index.get(other_index).copied())
                .collect();
            let found_index;

            if let Some(other_column_index) = other_columns.iter().min().copied() {
                found_index = other_column_index;
            } else {
                found_index = current_column;
                current_column += 1
            }

            column_index.insert(index, found_index);
        }

        // Step 2: given the column indexes, populate the strings
        //      corresponding to the sequences inserted in the graph
        let mut labels = vec![];
        let mut alignment_strings = vec![];
        let mut current_node_index_option: Option<NodeIndex>;

        for (label, &start) in self.labels.iter().zip(&self.start_indices) {
            labels.push(label.clone());
            current_node_index_option = Some(start);
            let mut item_list: Vec<String> =
                (0..current_column).map(|_| String::from("-")).collect();

            while let Some(current_node_index) = current_node_index_option {
                current_node_index_option = None;
                item_list[column_index[&current_node_index] as usize]
                    .clone_from(&self.graph[current_node_index].item);
                // iterate all out going edges
                for edge in self
                    .graph
                    .edges_directed(current_node_index, Direction::Outgoing)
                {
                    // found the edge with same label
                    if edge.weight().labels.iter().any(|item| item == label) {
                        current_node_index_option = Some(edge.target());
                        break;
                    }
                }
            }

            alignment_strings.push(item_list);
        }
        labels.into_iter().zip(alignment_strings).collect()
    }
}

#[cfg(test)]
mod consensus_tests {
    use crate::{alignment::SeqGraphAlignment, graph::POAGraph};

    use super::Consensus;

    #[test]
    fn consensus_test() {
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
        let seq3 = vec!["T".to_string(), "G".to_string(), "X".to_string()];
        let mut graph = POAGraph::new("seq_1".to_string(), seq1);
        let sg_aln = SeqGraphAlignment::align_seq_to_graph("seq_2".to_string(), seq2, &graph.graph);
        graph.add_alignment(sg_aln);
        let sg_aln = SeqGraphAlignment::align_seq_to_graph("seq_3".to_string(), seq3, &graph.graph);
        graph.add_alignment(sg_aln);
        let con = Consensus::new(graph.graph, graph.start_indices, graph.labels);
        let alns = con.compute();
        alns.iter().for_each(|v| {
            println!("{:?}", v);
        });
    }
}
