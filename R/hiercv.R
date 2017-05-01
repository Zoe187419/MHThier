#' Critical Function for Hierarchical FDR Controlling Procedure under Positive Dependence
#'
#' Given some parameters of the descendant hypotheses, return the critical function using the generalized hierarchical FDR controlling procedure under positive dependence (see Theorem 1 in Lynch et al. (2016)).
#'
#'@usage
#'  PositiveDeptCF(isTested, mi, li, l, alpha, rOffset)
#'@param isTested logical; if  \code{TRUE} (default), then the i-th hypotheses \eqn{H_i} will be tested; otherwise, \eqn{H_i} will not be tested.
#'@param mi the cardinality of the set of descendant hypotheses \eqn{H_i}.
#'@param li the number of leaf hypotheses in the set of descendant hypotheses \eqn{H_i}.
#'@param l the total number of leaf hypotheses.
#'@param alpha the significant level used to calculate the critical values to make decisions.
#'@param rOffset the offset increment for the number of rejections.
#'@return
#' A critical function of the index \eqn{i} and rejection number \eqn{R}.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W. (2016).
#'  On Procedures Controlling the FDR for Testing Hierarchically Ordered Hypotheses.
#'  \emph{arXiv preprint} arXiv:1612.04467.
#'@export
PositiveDeptCF <- function (isTested, mi, li, l, alpha, rOffset) {
  function (i, R) {
    if (isTested[i]) {
      li[i]*alpha/l * (mi[i]+R+rOffset-1)/mi[i]
    } else {
      0
    }
  }
}

#' Critical Function for Hierarchical FDR Controlling Procedure under Arbitrary Dependence.
#'
#' Given some parameters of the descendant hypotheses and ancestor hypotheses, return the critical function using the generalized hierarchical FDR controlling procedure under arbitrary dependence (see Theorem 2 in Lynch et al. (2016)).
#'
#'@usage
#'  ArbitraryDeptCF(isTested, mi, li, l, depth, gdi, alpha, rOffset)
#'@param isTested logical; if  \code{TRUE} (default), then the i-th hypotheses \eqn{H_i} will be tested; otherwise, \eqn{H_i} will not be tested.
#'@param mi the cardinality of the set of descendant hypotheses of \eqn{H_i}.
#'@param li the number of leaf hypotheses in the set of descendant hypotheses \eqn{H_i}.
#'@param l the total number of leaf hypotheses.
#'@param depth the cardinality of the set of ancestor hypotheses of \eqn{H_i}, referred to as the depth of \eqn{H_i}.
#'@param gdi the cardinality of the union set of all hypotheses with all depth no deeper than the given depth.
#'@param alpha the significant level used to calculate the critical values to make decisions.
#'@param rOffset the offset increment for the number of rejections.
#'@return
#' A critical function of the index \eqn{i} and rejection number \eqn{R}.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W. (2016).
#'  On Procedures Controlling the FDR for Testing Hierarchically Ordered Hypotheses.
#'  \emph{arXiv preprint} arXiv:1612.04467.
#'@export
ArbitraryDeptCF <- function (isTested, mi, li, l, depth, gdi, alpha, rOffset) {
    function (i, R) {
      if (isTested[i]) {
        li[i]*alpha/l * (mi[i]+R+rOffset-1)/mi[i] * 1/(sum(1 / (depth:(gdi-1) + mi[i]))+1)
      } else {
        0
      }
    }
  }

#' Critical Function for Hierarchical FDR Controlling Procedure under Block Positive Dependence
#'
#' Given some parameters of the descendant hypotheses, return the critical function using the generalized hierarchical FDR controlling procedure under block positive dependence (see Theorem 3 in Lynch et al. (2016)).
#'
#'@usage
#'  BlockPositiveCF(isTested, mi, li, l, alpha, rOffset)
#'@param isTested logical; if  \code{TRUE} (default), then the i-th hypotheses \eqn{H_i} will be tested; otherwise, \eqn{H_i} will not be tested.
#'@param mi the cardinality of the set of descendant hypotheses \eqn{H_i}.
#'@param li the number of leaf hypotheses in the set of descendant hypotheses \eqn{H_i}.
#'@param l the total number of leaf hypotheses.
#'@param alpha the significant level used to calculate the critical values to make decisions.
#'@param rOffset the offset increment for the number of rejections.
#'@return
#' A critical function of the index \eqn{i} and rejection number \eqn{R}.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W. (2016).
#'  On Procedures Controlling the FDR for Testing Hierarchically Ordered Hypotheses.
#'  \emph{arXiv preprint} arXiv:1612.04467.
#'@export
BlockPositiveCF <- function (isTested, mi, li, l, alpha, rOffset) {
 function (i, R) {
	if (isTested [i]) {
		if (mi[i] == 1) {
			(R+rOffset)*alpha / l
		} else {
			li[i]*(R+rOffset)*alpha / (l + li[i]*(R+rOffset-1)*alpha)
		}
	} else {
		0
	}
 }
}

#' Critical Function for Hierarchical FDR Controlling Procedure under Block Arbitrary Dependence.
#'
#' Given some parameters of the descendant hypotheses and ancestor hypotheses, return the critical function using the generalized hierarchical FDR controlling procedure under block arbitrary dependence (see Theorem 4 in Lynch et al. (2016)).
#'
#'@usage
#'  BlockArbitraryCF(isTested, mi, li, l, depth, fdi, alpha, rOffset)
#'@param isTested logical; if  \code{TRUE} (default), then the i-th hypotheses \eqn{H_i} will be tested; otherwise, \eqn{H_i} will not be tested.
#'@param mi the cardinality of the set of descendant hypotheses of \eqn{H_i}.
#'@param li the number of leaf hypotheses in the set of descendant hypotheses \eqn{H_i}.
#'@param l the total number of leaf hypotheses.
#'@param depth the cardinality of the set of ancestor hypotheses of \eqn{H_i}, referred to as the depth of \eqn{H_i}.
#'@param fdi the cardinality of the set of all hypotheses with the given depth.
#'@param alpha the significant level used to calculate the critical values to make decisions.
#'@param rOffset the offset increment for the number of rejections.
#'@return
#' A critical function of the index \eqn{i} and rejection number \eqn{R}.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W. (2016).
#'  On Procedures Controlling the FDR for Testing Hierarchically Ordered Hypotheses.
#'  \emph{arXiv preprint} arXiv:1612.04467.
#'@export
BlockArbitraryCF <- function (isTested, mi, li, l, depth, fdi, alpha, rOffset) {
    function (i, R) {
      if (isTested [i]) {
        if (mi[i] == 1) {
          (R+rOffset)*alpha / l * 1/(1 + sum(1/((1+depth):(fdi+depth-1))))
        } else {
          c = sum(sapply(1:(fdi-1), function(r) { (l - li[i]*alpha) / ((r+depth) * (l + li[i] * (r+depth-2)*alpha)) } ));
          li[i]*(R+rOffset)*alpha / (l + li[i]*(R+rOffset-1)*alpha) * 1 / (1+c)
        }
      } else {
        0
      }
    }
  }

#' Critical Function for all Hierarchical FDR Controlling Procedure under Various Types of Dependence.
#'
#' Given some parameters of the descendant hypotheses and ancestor hypotheses, return the critical function using the generalized hierarchical FDR controlling procedure under various type of dependence.
#'
#'@usage
#'  hierFDR.CF(isTested, mi, li, l, depth, fdi, gdi, alpha, rOffset, type)
#'@param isTested logical; if  \code{TRUE} (default), then the i-th hypotheses \eqn{H_i} will be tested; otherwise, \eqn{H_i} will not be tested.
#'@param mi the cardinality of the set of descendant hypotheses of \eqn{H_i}.
#'@param li the number of leaf hypotheses in the set of descendant hypotheses \eqn{H_i}.
#'@param l the total number of leaf hypotheses.
#'@param depth the cardinality of the set of ancestor hypotheses of \eqn{H_i}, referred to as the depth of \eqn{H_i}.
#'@param fdi the cardinality of the set of all hypotheses with the given depth.
#'@param gdi the cardinality of the union set of all hypotheses with all depth no deeper than the given depth.
#'@param alpha the significant level used to calculate the critical values to make decisions.
#'@param rOffset the offset increment for the number of rejections.
#'@param type the type of dependence structure of the hierarchically ordered hypotheses. Currently, we provide four types of dependence: \code{"positive"}, \code{"arbitrary"}, \code{"block positive"} and \code{"block arbitrary"}.
#'@return
#' A critical function of the index \eqn{i} and rejection number \eqn{R}.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W. (2016).
#'  On Procedures Controlling the FDR for Testing Hierarchically Ordered Hypotheses.
#'  \emph{arXiv preprint} arXiv:1612.04467.
#'@seealso \code{\link{PositiveDeptCF}}, \code{\link{ArbitraryDeptCF}}, \code{\link{BlockPositiveCF}}, \code{\link{BlockArbitraryCF}}
#'@export
hierFDR.CF <- function (isTested, mi, li, l, depth, fdi, gdi, alpha, rOffset, type = "block positive") {
  if (type == "positive") {
    PositiveDeptCF (isTested, mi, li, l, alpha, rOffset);
  } else if (type == "arbitrary") {
    ArbitraryDeptCF (isTested, mi, li, l, depth, gdi, alpha, rOffset);
  } else if (type == "block positive") {
    BlockPositiveCF (isTested, mi, li, l, alpha, rOffset);
  } else if (type == "block arbitrary") {
    BlockArbitraryCF (isTested, mi, li, l, depth, fdi, alpha, rOffset);
  }
}


#' Decision-making Function for Generalized Stepwise Procedures
#'
#' Given a set of p-values, return the decisions using the generalized stepwise procedure.
#'
#'@usage
#'  stepwise(p, CF, k)
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param CF critical functions of the corresponding multiple testing procedure
#'@param k the order of a generalized stepwise procedure.
#'@return
#' logical values of each hypothesis being rejected or not, if \code{TRUE}, then the hypothesis is rejected; otherwise, the hypothesis is not rejected.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W. (2016).
#'  On Procedures Controlling the FDR for Testing Hierarchically Ordered Hypotheses.
#'  \emph{arXiv preprint} arXiv:1612.04467.
#'
#'  Tamhane, A. C., Liu, W., & Dunnett, C. W. (1998).
#'  A generalized step-up-down multiple test procedure.
#'  \emph{The Canadian Journal of Statistics/La Revue Canadienne de Statistique}, 353-363.
#'
#'  Sarkar, S. K. (2002).
#'  Some results on false discovery rate in stepwise multiple testing procedures.
#'  \emph{Annals of statistics}, 239-257.
#'
#'@export
stepwise <- function (p, CF, k) {
  tot = length(p);
  r = k;
  crits = sapply(1:tot, function (i) { CF(i,r) });
  s = sum(p <= crits);
  doStepA = (r > s);
  for (i in (1:tot)) {
    if (doStepA) {
      r = s;
    } else {
      r = s+1;
    }
    crits = sapply(1:tot, function (i) { CF(i,r) });
    s = sum(p <= crits);
    if ((doStepA) && (r <= s)) {
      break;
    } else if ((!doStepA) && (r > s)) {
      r = r-1;
      crits = sapply(1:tot, function (i) { CF(i,r) });
      break;
    }
  }
  p <= crits
}

#' Decision-making Function for Generalized Hierarchical Testing Procedures
#'
#' Given a set of p-values, return the decisions using the generalized stepwise procedure.
#'
#'@usage
#'  hier.test(tree, pvals, alpha, type)
#'@import igraph
#'@param tree the edgelist parameterizing the hierarchical structure between hypotheses. The edges must be stored so that each edge is a row of a two column matrix, where the first column gives the parent and the second gives the child.
#'@param pvals a vector of raw p-values resulting from an experiment. The names of this vector should be contained in the edgelist parameterizing the hierarchical structure between hypothesis, inputted as \code{tree}
#'@param alpha the significant level used to calculate the critical values to make decisions.
#'@param type the type of dependence structure of the hierarchically ordered hypotheses. Currently, we provide four types of dependence: \code{"positive"}, \code{"arbitrary"}, \code{"block positive"} and \code{"block arbitrary"}.
#'@return
#' logical values of each hypothesis being rejected or not, if \code{TRUE}, then the hypothesis is rejected; otherwise, the hypothesis is not rejected.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W. (2016).
#'  On Procedures Controlling the FDR for Testing Hierarchically Ordered Hypotheses.
#'  \emph{arXiv preprint} arXiv:1612.04467.
#'
#'@examples
#' library(igraph)
#' library(ape)
#' library(structSSI)
#' library(phyloseq)
#' data("chlamydiae")
#' environments <- sample_data(chlamydiae)$SampleType
#' abundances <- otu_table(chlamydiae)
#'
#' graph.tree <- as.igraph(phy_tree(chlamydiae))
#' edge.tree <- get.edgelist(graph.tree)
#' pVal <- treePValues(edge.tree, abundances, environments)
#' pVal[which(is.na(pVal))] = 1;		# these have all 0 abundances in every environment
#' decision1 <- hier.test(tree = graph.tree, pvals = pVal, alpha = 0.01, type = "positive")
#' decision2 <- hier.test(tree = graph.tree, pvals = pVal, alpha = 0.01, type = "arbitrary")
#' ## show the number of rejections under different types of dependence
#' length(which(decision1 == TRUE)); length(which(decision2 == TRUE))
#'@export
hier.test <- function (tree, pvals, alpha, type = "block positive") {

  # get the hypotheses which do not appear in the second column of edgelist
  m = length(V(tree));
  edgelist <- get.edgelist(tree, names = FALSE)
  topHyps = (1:m)[-edgelist[,2]];

  # can't use shortest.paths since for large trees it runs out of memory
  V(tree)$depth = shortest.paths(tree, to=topHyps)+1;	#assign the depth to each hypothesis

  #graphDiameter <- diameter(tree);
  #descendants = topHyps;
  #V(tree)[descendants]$depth = 1;
  #for (d in (2:graphDiameter)) {
  #	newDescendants <- unlist(neighborhood(tree, nodes = descendants, mode = "out", order = 1));
  #	descendants = newDescendants[-match(descendants, newDescendants)]
  #	V(tree)[descendants]$depth = d;
  #}

  # determine the metadata such as total hypotheses in the subtree, leafs, etc
  for(d in (max(V(tree)$depth):1)) {
    inds = which(V(tree)$depth == d)
    for (i in (1:length(inds))) {
      childInds = edgelist[edgelist[,1] == inds[i],2];
      if (length(childInds) == 0) {
        V(tree)[inds[i]]$li = 1;
        V(tree)[inds[i]]$mi = 1;
      } else {
        V(tree)[inds[i]]$li = sum(V(tree)[childInds]$li);
        V(tree)[inds[i]]$mi = 1 + sum(V(tree)[childInds]$mi);
      }
    }
  }
  l = sum(V(tree)[topHyps]$li);

  rOffset = 0;
  isTested = rep(FALSE, m);

  # test each family 1 through d
  for(d in (1:max(V(tree)$depth))) {
    inds = which(V(tree)$depth == d);		# get hypotheses in family d
    gdi = length(which(V(tree)$depth <= d));
    if (d == 1) {
      # test everything in the first family
      isTested [inds] = rep(TRUE, length(inds));
    } else {
      # only test the hypotheses whose parent was rejected
      parentInds = sapply(inds, function(x) edgelist[edgelist[,2] == x,1]);
      isTested [inds] = V(tree)[parentInds]$rejected;
    }
    crits = hierFDR.CF(isTested [inds], V(tree)[inds]$mi, V(tree)[inds]$li, l, d, length(inds), gdi, alpha, rOffset, type)

    # run the empirical stepwise procedure initialized to the maximum
    V(tree)[inds]$rejected = stepwise(pvals[inds], crits, length(inds));
    rOffset = rOffset + sum(V(tree)[inds]$rejected);
  }

  # these are accepted by default but more specifically not tested
  V(tree)[! isTested]$rejected = NA;
  V(tree)$rejected
}

