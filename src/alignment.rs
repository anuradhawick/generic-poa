use ndarray::Array2;
use petgraph::{
    algo::toposort,
    graph::{DiGraph, NodeIndex},
    visit::EdgeRef,
    Direction,
};
use std::cmp::max;
use std::collections::HashMap;

use crate::graph::{EdgeData, NodeData};

const MATCHSCORE: i32 = 1;
const MISMATCHSCORE: i32 = -1;
const GAP: i32 = -2;

type Matrix = Array2<i32>;

#[derive(Debug)]
enum Move {
    Insertion,
    Deletion,
    Match,
}

#[derive(Debug, Eq)]
struct Candidate {
    score: i32,
    graph_pos: usize,
    seq_pos: usize,
}

// Implementing PartialEq is necessary to use PartialOrd
impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool {
        (self.score, self.graph_pos, self.seq_pos) == (other.score, other.graph_pos, other.seq_pos)
    }
}

// Implementing Ord to define complete comparison
impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.score, self.graph_pos, self.seq_pos).cmp(&(
            other.score,
            other.graph_pos,
            other.seq_pos,
        ))
    }
}

impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

pub struct SeqGraphAlignment {
    pub seq: Vec<String>,
    // pub graph: DiGraph<NodeData, EdgeData>,
    pub seq_match_positions: Vec<Option<i32>>,
    pub graph_match_node_indices: Vec<Option<NodeIndex>>,
}

impl SeqGraphAlignment {
    pub fn align_seq_to_graph(seq: Vec<String>, graph: &DiGraph<NodeData, EdgeData>) -> Self {
        let indices = toposort(&graph, None).unwrap();
        let (
            mut scores,
            mut backtrack_scores_seq,
            mut backtrack_score_graph,
            node_index_to_matrix_pos,
            matrix_pos_to_node_index,
        ) = SeqGraphAlignment::prep_dp_matrix(&graph, &indices, &seq);

        // println!("Starting with scores \n {:?}", scores);

        // alignment
        for (i, &index) in indices.iter().enumerate() {
            let pbase = graph[index].item.as_str();

            // println!("i={i:?}");

            for (j, sbase) in seq.iter().enumerate() {
                let mut candidates = vec![Candidate {
                    score: scores[[i + 1, j]] + GAP,
                    graph_pos: i + 1,
                    seq_pos: j,
                }];
                // candidates = [(scores[i+1, j] + self._gap, i+1, j, "INS")]
                let mut prev_positions: Vec<i32> = graph
                    .edges_directed(index, Direction::Incoming)
                    .map(|e| node_index_to_matrix_pos[&e.source()] as i32)
                    .collect();

                if prev_positions.is_empty() {
                    prev_positions = vec![-1];
                }

                for prev_node_pos in prev_positions {
                    // candidates += [(scores[predIndex+1, j+1] + self._gap, predIndex+1, j+1, "DEL")]
                    candidates.push(Candidate {
                        score: scores[[(prev_node_pos + 1) as usize, j + 1]] + GAP,
                        graph_pos: (prev_node_pos + 1) as usize,
                        seq_pos: j + 1,
                    });
                    // candidates += [(scores[predIndex+1, j] + self.matchscore(sbase, pbase), predIndex+1, j, "MATCH")]
                    candidates.push(Candidate {
                        score: scores[[(prev_node_pos + 1) as usize, j]]
                            + if sbase == pbase {
                                MATCHSCORE
                            } else {
                                MISMATCHSCORE
                            },
                        graph_pos: (prev_node_pos + 1) as usize,
                        seq_pos: j,
                    });
                }

                // println!("\tj={}", j);
                // println!("\tcandidates={:?}", candidates);
                let chosen_candidate = candidates.iter().max().unwrap();

                // println!("\t\tchosen candidate={:?}", chosen_candidate);
                scores[[i + 1, j + 1]] = chosen_candidate.score;
                backtrack_score_graph[[i + 1, j + 1]] = chosen_candidate.graph_pos as i32;
                backtrack_scores_seq[[i + 1, j + 1]] = chosen_candidate.seq_pos as i32;
            }
        }

        // println!("scores \n{:?}", scores);
        // println!("backtrack_score_graph \n{:?}", backtrack_score_graph);
        // println!("backtrack_scores_seq \n{:?}", backtrack_scores_seq);

        let (seq_match_positions, graph_match_node_indices) = SeqGraphAlignment::backtrack(
            &graph,
            &indices,
            &scores,
            &backtrack_scores_seq,
            &backtrack_score_graph,
            &matrix_pos_to_node_index,
        );

        Self {
            seq,
            // graph,
            seq_match_positions,
            graph_match_node_indices,
        }
    }

    pub fn to_string(&self, graph: &DiGraph<NodeData, EdgeData>) -> (String, String, String) {
        let s1: Vec<String> = self
            .seq_match_positions
            .iter()
            .map(|item| {
                if let Some(pos) = item {
                    self.seq[*pos as usize].to_string()
                } else {
                    String::from("-")
                }
            })
            .collect();

        let s2: Vec<String> = self
            .graph_match_node_indices
            .iter()
            .map(|item| {
                if let Some(pos) = item {
                    graph[*pos].item.clone()
                } else {
                    String::from("-")
                }
            })
            .collect();

        let m: Vec<String> = s1
            .iter()
            .zip(&s2)
            .map(|(f1, f2)| {
                if f1 == f2 {
                    String::from("|")
                } else {
                    String::from(" ")
                }
            })
            .collect();

        (s1.join(""), m.join(""), s2.join(""))
    }

    fn backtrack(
        graph: &DiGraph<NodeData, EdgeData>,
        indices: &[NodeIndex],
        scores: &Matrix,
        backtrack_scores_seq: &Matrix,
        backtrack_score_graph: &Matrix,
        matrix_pos_to_node_index: &HashMap<usize, NodeIndex>,
    ) -> (Vec<Option<i32>>, Vec<Option<NodeIndex>>) {
        let shape = scores.shape();
        let mut besti = shape[0] as i32 - 1;
        let mut bestj = shape[1] as i32 - 1;

        // global alignment scenario
        if true {
            let terminal_indices: Vec<i32> = indices
                .iter()
                .enumerate()
                .filter_map(|(pos, index)| {
                    if graph.edges_directed(*index, Direction::Outgoing).count() == 0 {
                        Some(pos as i32)
                    } else {
                        None
                    }
                })
                .collect();
            besti = terminal_indices[0] + 1;
            let mut best_score = scores[[besti as usize, bestj as usize]];

            for &i in &terminal_indices[1..] {
                let score = scores[[(i + 1) as usize, bestj as usize]];
                if score > best_score {
                    best_score = score;
                    besti = i + 1;
                }
            }
            println!("Terminal indices= {:?}", terminal_indices);
            println!("besti= {:?}", besti);
            println!("bestj= {:?}", bestj);
            println!("best score= {:?}", best_score);
        }

        let mut graph_match_node_indices = vec![];
        let mut seq_match_positions = vec![];
        let mut nexti;
        let mut nextj;
        let mut current_seq_pos;
        let mut current_node_index;
        let is_global = true;

        while (is_global || scores[[besti as usize, bestj as usize]] > 0)
            && (besti != 0 || bestj != 0)
        {
            nexti = backtrack_score_graph[[besti as usize, bestj as usize]];
            nextj = backtrack_scores_seq[[besti as usize, bestj as usize]];

            println!("-----");
            println!("nexti = {nexti}");
            println!("nextj = {nextj}");
            current_seq_pos = bestj - 1;
            current_node_index = if besti > 0 {
                matrix_pos_to_node_index[&((besti - 1) as usize)]
            } else {
                NodeIndex::from(u32::MAX)
            };

            seq_match_positions.insert(
                0,
                if nextj != bestj {
                    Some(current_seq_pos)
                } else {
                    None
                },
            );
            graph_match_node_indices.insert(
                0,
                if nexti != besti {
                    Some(current_node_index)
                } else {
                    None
                },
            );
            besti = nexti;
            bestj = nextj;
        }

        println!("{:?}", seq_match_positions);
        println!("{:?}", graph_match_node_indices);

        (seq_match_positions, graph_match_node_indices)
    }

    fn prep_dp_matrix(
        graph: &DiGraph<NodeData, EdgeData>,
        indices: &[NodeIndex],
        seq: &[String],
    ) -> (
        Matrix,
        Matrix,
        Matrix,
        HashMap<NodeIndex, usize>,
        HashMap<usize, NodeIndex>,
    ) {
        let l1 = graph.node_count();
        let l2 = seq.len();
        let mut scores = Array2::<i32>::zeros((l1 + 1, l2 + 1));
        let mut node_index_to_matrix_pos: HashMap<NodeIndex, usize> = HashMap::new();
        let mut matrix_pos_to_node_index: HashMap<usize, NodeIndex> = HashMap::new();
        let backtrack_scores_seq = Array2::<i32>::zeros((l1 + 1, l2 + 1));
        let backtrack_score_graph = Array2::<i32>::zeros((l1 + 1, l2 + 1));

        for (position, &node_index) in indices.iter().enumerate() {
            node_index_to_matrix_pos.insert(node_index, position);
            matrix_pos_to_node_index.insert(position, node_index);
        }

        // if global alignment
        if true {
            for i in 0..l2 + 1 {
                scores[[0, i]] = (i as i32) * GAP; // mul by gap later
            }

            for (matrix_pos, &node_index) in indices.iter().enumerate() {
                let prev_node_indices = graph.edges_directed(node_index, Direction::Incoming);
                let mut best_score = 0;
                for edge in prev_node_indices {
                    // if score is zero and there are incoming edges (we're in this loop), it is not possible
                    if best_score == 0 {
                        best_score = scores[[node_index_to_matrix_pos[&edge.source()] + 1, 0]];
                    }
                    best_score = max(
                        best_score,
                        scores[[node_index_to_matrix_pos[&edge.source()] + 1, 0]],
                    );
                }
                scores[[matrix_pos + 1, 0]] = best_score + GAP;
            }
        }
        (
            scores,
            backtrack_scores_seq,
            backtrack_score_graph,
            node_index_to_matrix_pos,
            matrix_pos_to_node_index,
        )
    }
}

#[cfg(test)]
mod graph_tests {
    use crate::graph::POAGraph;

    use super::SeqGraphAlignment;

    #[test]
    fn align_seq_test() {
        let seq1 = vec![
            // "A".to_string(),
            "M".to_string(),
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
        let graph = POAGraph::new(seq1);
        let sg_aln = SeqGraphAlignment::align_seq_to_graph(seq2, &graph.graph);
        let (f1, m, f2) = sg_aln.to_string(&graph.graph);
        println!("{f1}");
        println!("{m}");
        println!("{f2}");
    }
}
