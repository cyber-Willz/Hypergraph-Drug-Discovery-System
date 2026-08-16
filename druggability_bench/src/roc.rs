//! ROC-AUC via the rank-sum (Mann-Whitney U) identity, which is exact and
//! avoids numerically-integrating a thresholded ROC curve:
//!
//! ```text
//! AUC = (sum_of_ranks(positive_scores) - n_pos*(n_pos+1)/2) / (n_pos * n_neg)
//! ```
//!
//! where ranks are 1-indexed over all `n_pos + n_neg` scores pooled
//! together, ascending, with tied scores receiving the *average* rank of
//! their tied block (the standard tie-correction; without it, AUC is
//! biased on data with many repeated scores, which classical centrality
//! baselines like degree produce a lot of on small contact graphs).
//!
//! This is equivalent to: AUC = P(score(random positive) > score(random
//! negative)) + 0.5 * P(tie), which is the quantity the ROC curve's area
//! actually represents.

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AucResult {
    pub auc: f64,
    pub n_pos: usize,
    pub n_neg: usize,
}

#[derive(Debug, thiserror::Error, PartialEq)]
pub enum RocError {
    #[error("need at least one positive and one negative label to define an ROC curve (got {n_pos} positive, {n_neg} negative)")]
    DegenerateLabels { n_pos: usize, n_neg: usize },
    #[error("scores and labels must have the same length (got {n_scores} scores, {n_labels} labels)")]
    LengthMismatch { n_scores: usize, n_labels: usize },
}

/// Compute ROC-AUC for `scores` against binary ground-truth `labels`
/// (`true` = positive class, e.g. "this residue lines a known cryptic
/// pocket"). Higher `scores` should indicate more likely positive.
pub fn roc_auc(scores: &[f64], labels: &[bool]) -> Result<AucResult, RocError> {
    if scores.len() != labels.len() {
        return Err(RocError::LengthMismatch { n_scores: scores.len(), n_labels: labels.len() });
    }
    let n_pos = labels.iter().filter(|&&l| l).count();
    let n_neg = labels.len() - n_pos;
    if n_pos == 0 || n_neg == 0 {
        return Err(RocError::DegenerateLabels { n_pos, n_neg });
    }

    // Sort indices by score ascending, then assign average rank within
    // each block of tied scores.
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| scores[a].partial_cmp(&scores[b]).unwrap());

    let mut ranks = vec![0.0f64; scores.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i;
        while j + 1 < order.len() && scores[order[j + 1]] == scores[order[i]] {
            j += 1;
        }
        // Ranks are 1-indexed; average rank of the tied block [i, j].
        let avg_rank = ((i + 1) + (j + 1)) as f64 / 2.0;
        for &idx in &order[i..=j] {
            ranks[idx] = avg_rank;
        }
        i = j + 1;
    }

    let rank_sum_pos: f64 =
        labels.iter().zip(ranks.iter()).filter(|(&l, _)| l).map(|(_, &r)| r).sum();

    let u = rank_sum_pos - (n_pos * (n_pos + 1)) as f64 / 2.0;
    let auc = u / (n_pos * n_neg) as f64;
    Ok(AucResult { auc, n_pos, n_neg })
}

/// Pool scores and labels from multiple structures into one AUC (as
/// opposed to averaging per-structure AUCs). Standard in benchmarks with
/// small per-structure residue counts, where per-structure AUC is noisy.
pub fn pooled_roc_auc(runs: &[(Vec<f64>, Vec<bool>)]) -> Result<AucResult, RocError> {
    let mut scores = Vec::new();
    let mut labels = Vec::new();
    for (s, l) in runs {
        scores.extend_from_slice(s);
        labels.extend_from_slice(l);
    }
    roc_auc(&scores, &labels)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn perfect_separation_gives_auc_one() {
        let scores = vec![0.1, 0.2, 0.8, 0.9];
        let labels = vec![false, false, true, true];
        let r = roc_auc(&scores, &labels).unwrap();
        assert_relative_eq!(r.auc, 1.0);
    }

    #[test]
    fn inverted_separation_gives_auc_zero() {
        let scores = vec![0.9, 0.8, 0.2, 0.1];
        let labels = vec![false, false, true, true];
        let r = roc_auc(&scores, &labels).unwrap();
        assert_relative_eq!(r.auc, 0.0);
    }

    #[test]
    fn all_tied_scores_give_auc_half() {
        let scores = vec![0.5, 0.5, 0.5, 0.5];
        let labels = vec![true, false, true, false];
        let r = roc_auc(&scores, &labels).unwrap();
        assert_relative_eq!(r.auc, 0.5);
    }

    #[test]
    fn matches_brute_force_pair_counting() {
        // Cross-check the rank-sum formula against the O(n_pos*n_neg)
        // definition directly, on a case with partial ties.
        let scores = vec![1.0, 2.0, 2.0, 3.0, 1.5, 4.0];
        let labels = vec![false, true, false, true, true, false];
        let r = roc_auc(&scores, &labels).unwrap();

        let pos: Vec<f64> =
            scores.iter().zip(&labels).filter(|(_, &l)| l).map(|(&s, _)| s).collect();
        let neg: Vec<f64> =
            scores.iter().zip(&labels).filter(|(_, &l)| !l).map(|(&s, _)| s).collect();
        let mut total = 0.0;
        for &p in &pos {
            for &n in &neg {
                if p > n {
                    total += 1.0;
                } else if p == n {
                    total += 0.5;
                }
            }
        }
        let brute_force_auc = total / (pos.len() * neg.len()) as f64;
        assert_relative_eq!(r.auc, brute_force_auc);
    }

    #[test]
    fn degenerate_labels_rejected() {
        let scores = vec![0.1, 0.2, 0.3];
        let labels = vec![true, true, true];
        assert_eq!(
            roc_auc(&scores, &labels),
            Err(RocError::DegenerateLabels { n_pos: 3, n_neg: 0 })
        );
    }
}
