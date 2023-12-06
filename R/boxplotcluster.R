#' Boxplot Cluster Function
#'
#' @description The function \code{boxplotcluster} implements a special clustering method
#' based on boxplot statistics. Following Arroyo-Maté-Roque (2006), the function calculates
#' the distance between rows or columns of the dataset using the generalized Minkowski metric as described
#' by Ichino and Yaguchi (1994). The distance measure gives more weight to differences between
#' quartiles than to differences between extremes, making it less sensitive to outliers. Further,
#' the function calculates the silhouette width for different numbers of clusters (Rousseeuw 1987)
#' and selects the number of clusters that maximizes the average silhouette width
#' (unless a specific number of clusters is provided by the user).\cr
#' Visit this \href{https://drive.google.com/file/d/1Qb9r8amiSa1OMmCx3ngpeXAGT7kebInq/view?usp=share_link}{LINK} to access the package's vignette.\cr
#'
#' @param x A dataframe representing the input dataset (see below and the \code{Details} section).
#' @param target.var A character vector specifying the name of the numerical variable in the input dataset that
#' will be used for analysis. If the input dataset is in \code{wide} format, this can be left as default (\code{NULL}). If the dataset is
#' in \code{long} format, this parameter has to indicate the variable containing the data values that will
#' be grouped by the different levels of the variable fed via the \code{group.var} parameter.
#' @param group.var A character vector specifying the grouping variable. It can be left as default (\code{NULL}) if the
#' input dataset is in \code{wide} format.
#' @param calc.type A string specifying the units to be clustered if the input dataset is in wide format (either \code{columns} or \code{rows}; the former is the default).
#' @param aggl.meth A string specifying the agglomeration method to be used in hierarchical
#' clustering. Defaults to "ward.D2". For other methods see \code{\link[stats]{hclust}}.
#' @param part An optional integer specifying the desired number of clusters. If not provided,
#' the function selects the number of clusters that maximises the average silhouette width.
#' @param silh.col A logical value, which takes \code{TRUE} (default) or \code{FALSE} if the user wants to give colour to the
#' silhouette plot reflecting the cluster partition.
#' @param cex.dndr.lab A numeric specifying the character expansion factor for the labels in the
#' dendrogram plot. Defaults to 0.75.
#' @param cex.sil.lab A numeric specifying the character expansion factor for the labels in the
#' silhouette plot. Defaults to 0.75.
#' @param oneplot A logical value, which takes \code{TRUE} (default) or \code{FALSE} if the user wants or does not want the plots to be visualised
#' in a single window.
#'
#' @details The function first calculates the pairwise distance between each unit of the input dataset
#'using the Ichino-Yaguchi dissimilarity measure (equations 7 and 8 in Arroyo-Maté-Roque (2006)).
#' The distance between A and B is defined as:\cr
#'
#' \eqn{(0.5 * (abs(m1 - m2) + 2 * abs(q1 - q2) + 2 * abs(Me1 - Me2) + 2 * abs(Q1 - Q2) + abs(M1 - M2))) /4}\cr
#'
#'   where\cr
#'
#'   m1 <- min(A)\cr
#'   m2 <- min(B)\cr
#'   q1 <- quantile(A, probs = 0.25)\cr
#'   q2 <- quantile(B, probs = 0.25)\cr
#'   Q1 <- quantile(A, probs = 0.75)\cr
#'   Q2 <- quantile(B, probs = 0.75)\cr
#'   M1 <- max(A)\cr
#'   M2 <- max(B)\cr
#'   Me1 <- median(A)\cr
#'   Me2 <- median(B)\cr
#'
#' The distance matrix is then used to perform a hierarchical clustering. Also, the function calculates the
#' silhouette width for different numbers of clusters and selects the number of clusters
#' that maximises the average silhouette width (unless a specific number of clusters is provided by the user).\cr
#'
#' The silhouette method allows to measure how 'good' is the selected cluster solution. If the parameter \code{part}
#' is left empty (default), an optimal cluster solution is obtained. The optimal partition is selected via an iterative
#' procedure which identifies at which cluster solution the highest average silhouette width is achieved.
#' If a user-defined partition is needed, the user can input the desired number of clusters using the parameter \code{part}.
#' In either case, an additional plot is returned besides the cluster dendrogram and the silhouette plot; it displays a
#' scatterplot in which the cluster solution (x-axis) is plotted against the average silhouette width (y-axis).
#' A black dot represents the partition selected either by the iterative procedure or by the user.\cr
#'
#' In summary, the function generates a series of plots to visualise the results:\cr
#'
#' (a) boxplots colored by cluster membership,\cr
#' (b) a dendrogram (where clusters are indicated by rectangles whose color is consistent with the color assigned to the
#' boxplot in the previous plot),\cr
#' (c) a silhouette plot (with optional by-cluster colours), and\cr
#' (d) a plot of the average silhouette width vs. the number of clusters.\cr
#'
#' The silhouette plot is obtained via the \code{silhouette()} function out from the \code{cluster} package.
#' For a detailed description of the silhouette plot, its rationale, and its interpretation, see Rousseeuw 1987.
#'
#' Input dataset format:\cr
#' When the dataset is in \code{wide} format, each row corresponds to a distinct unit, and each column corresponds
#' to a different variable. In this format, it is typically assumed that all units have the same number of observations.
#' The element representing the units (either rows or columns) will be clustered.\cr
#'
#' When the dataset is in \code{long} format, it consists of rows representing individual observations. One column indicates
#' the variable name (grouping variable) and another column contains the measurements values, or viceversa. If the input dataset
#' comes in this format, the units to be clustered are created by grouping the observations (i.e., rows of the
#' dataframe) by \code{group.var}. If the input dataset is in \code{long} format, groups can feature a different number of observations.\cr
#'
#' For actual examples of both formats, see the \code{Examples} section below.
#'
#'
#'@return The function returns a list storing the following components
#' \itemize{
##'  \item{\code{distance.matrix:}} distance matrix reporting the distance values.
##'  \item{\code{units.by.cluster:}} a list of the input dataset's units, grouped by cluster membership.
##'  \item{\code{avr.silh.width.by.n.of.clusters:}} average silhouette width by number of clusters.
##'  \item{\code{partition.silh.data:}} silhouette data for the selected partition.
##' }
#'
#' @references Arroyo, J., Maté, C., & Roque, A. M-S. (2006). Hierarchical Clustering for Boxplot Variables.
#' In Studies in Classification, Data Analysis, and Knowledge Organization (pp. 59–66).
#' Springer Berlin Heidelberg.
#'
#' @references Ichino, M., & Yaguchi, H. (1994). Generalized Minkowski Metrics for Mixed
#' Feature-Type Data Analysis. IEEE Transactions on Systems, Man, and Cybernetics, 24(4), 698-708.
#'
#' @references Rousseeuw, P J. (1987). Silhouettes: A graphical aid to the interpretation and validation of cluster analysis,
#' Journal of Computational and Applied Mathematics 20, 53-65.
#'
#'
#' @examples
#'
#' ## EXAMPLE 1
#'
#' # Create a toy dataset in WIDE format
#' df <- data.frame(
#' a = rnorm(30, mean = 30, sd = 5),
#' b = rnorm(30, mean = 40, sd = 5),
#' c = rnorm(30, mean = 20, sd = 5),
#' d = rnorm(30, mean = 25, sd = 5),
#' e = rnorm(30, mean = 50, sd = 5),
#' f = rnorm(30, mean = 10, sd = 5),
#' g = rnorm(30, mean = 100, sd = 5),
#' h = rnorm(30, mean = 20, sd = 5),
#' i = rnorm(30, mean = 40, sd = 5),
#' l = rnorm(30, mean = 35, sd = 5),
#' m = rnorm(30, mean = 35, sd = 5),
#' n = rnorm(30, mean = 70, sd = 5),
#' o = rnorm(30, mean = 20, sd = 5),
#' p = rnorm(30, mean = 70, sd = 5),
#' q = rnorm(30, mean = 90, sd = 5)
#' )
#'
#' # Run the function
#' result <- boxplotcluster(df)
#'
#' # Same as above, but selecting a 4-cluster solution
#' result <- boxplotcluster(df, part=4)
#'
#' # Same as above, but the rows are clustered
#' result <- boxplotcluster(df, calc.type="rows", part=4)
#'
#'
#'
#' ## EXAMPLE 2
#' # Create a toy dataset in WIDE format, representing archaeological stone
#' # flake length (cm) by raw material
#'
#' df <- data.frame(
#' Basalt = c(7.0, 7.0, 7.7, 8.2, 10.3, 10.3, 10.3, 10.8, 11.0, 13.0, 13.9, 14.6, 1.0),
#' Chert = c(2.9, 4.8, 5.3, 5.8, 5.8, 6.2, 6.5, 7.7, 7.7, 7.9, 8.9, 9.6, 2.0),
#' Obsidian = c(2.2, 2.4, 3.1, 4.3, 5.0, 5.5, 5.8, 6.0, 6.2, 7.2, 7.4, 7.7, 2.0),
#' Quartzite = c(5.5, 5.5, 7.0, 7.4, 7.7, 7.9, 8.6, 8.9, 9.4, 9.6, 10.6, 10.8, 1.0),
#' Granite = c(4.0, 4.5, 6.0, 6.8, 7.0, 7.8, 8.1, 8.4, 9.0, 10.0, 10.5, 11.0, 1.0),
#' Sandstone = c(3.0, 3.2, 4.0, 4.5, 4.9, 5.2, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 1.0),
#' Limestone = c(3.5, 4.0, 4.8, 5.2, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 1.5),
#' Slate = c(2.0, 2.4, 3.0, 3.6, 4.0, 4.4, 5.0, 5.6, 6.0, 6.4, 7.0, 7.6, 1.2)
#' )
#'
#' # Run the function to cluster the columns (default); cluster solution is
#' # selected by the iterative method (default)
#'
#' result <- boxplotcluster(df)
#'
#'
#'
#' ## EXAMPLE 3
#' #  Create a toy dataset in LONG format
#'
#' n_units <- 20
#' n_groups <- 10
#' measurements_per_group <- 4
#'
#' long_data <- data.frame(
#' SubjectID = rep(paste0("Unit", 1:n_units), each = n_groups * measurements_per_group),
#' Grouping_var = rep(rep(paste0("M", 1:n_groups), each = measurements_per_group), n_units),
#' Value = runif(n_units * n_groups * measurements_per_group))
#'
#' #  Run the analysis, specifying the target variable and the grouping variable,
#' # and selecting a 3-cluster solution
#' result <- boxplotcluster(long_data, target.var = "Value", group.var = "Grouping_var", part=3)
#'
#'
#' @importFrom grDevices rainbow
#' @importFrom graphics abline axis boxplot layout par
#' @importFrom stats as.dist cutree hclust median quantile rect.hclust
#'
#' @export
#'
#' @seealso \code{\link[cluster]{silhouette}}, \code{\link[stats]{hclust}}
#'
#'
boxplotcluster <- function(x, target.var = NULL, group.var = NULL, calc.type = "columns", aggl.meth = "ward.D2", part = NULL, silh.col = TRUE, cex.dndr.lab = 0.75, cex.sil.lab = 0.75, oneplot=TRUE) {

  # Save current par settings
  oldpar <- par(no.readonly = TRUE)
  # Ensure settings are restored when function exits
  on.exit(par(oldpar))

  if(oneplot==TRUE){
    m <- rbind(c(1,2), c(3,4))
    layout(m)
  }

  # Convert data to a list of groups
  if (!is.null(target.var) && !is.null(group.var)) {
    # Long format: Split the data by group.var
    if (!target.var %in% names(x) || !group.var %in% names(x)) {
      stop("target.var or group.var not found in the dataframe")
    }
    groups <- split(x[[target.var]], x[[group.var]])
  } else {
    # Original format: Each column or row is treated as a group
    if (calc.type == "rows") {
      x <- t(x)
    }
    groups <- as.list(as.data.frame(x))
  }

  # Define the function for distance calculation based on box-plot statistics
  calc_dist <- function(group1, group2) {
    m1 <- min(group1)
    m2 <- min(group2)
    q1 <- quantile(group1, probs = 0.25)
    q2 <- quantile(group2, probs = 0.25)
    Q1 <- quantile(group1, probs = 0.75)
    Q2 <- quantile(group2, probs = 0.75)
    M1 <- max(group1)
    M2 <- max(group2)
    Me1 <- median(group1)
    Me2 <- median(group2)

    (0.5 * (abs(m1 - m2) +
              2 * abs(q1 - q2) +
              2 * abs(Me1 - Me2) +
              2 * abs(Q1 - Q2) +
              abs(M1 - M2))) / 4
  }

  # Calculate the pairwise distances between groups
  num_groups <- length(groups)
  d <- matrix(NA, nrow = num_groups, ncol = num_groups)
  group_labels <- names(groups)
  rownames(d) <- colnames(d) <- group_labels

  for (i in 1:num_groups) {
    for (j in 1:num_groups) {
      if (i != j) {
        d[i, j] <- calc_dist(groups[[i]], groups[[j]])
      } else {
        d[i, j] <- 0
      }
    }
  }

  # Convert the distance matrix to a proper labeled matrix
  # d<- as.matrix(d, labels = TRUE)
  # colnames(d) <- rownames(d) <- colnames(x)

  # Convert the distance matrix to a dist object
  d <- as.dist(d)

  # Perform hierarchical clustering
  fit <- hclust(d, method = aggl.meth)

  # Determine the maximum number of clusters
  max.ncl <- num_groups - 1

  # Initialize a vector to store the silhouette width values
  sil.width.val <- numeric(max.ncl - 1)

  # Initialize a vector for the possible number of clusters
  sil.width.step <- c(2:max.ncl)

  # Calculate the silhouette width for each possible number of clusters
  for (i in 2:max.ncl) {
    counter <- i - 1
    clust <- cluster::silhouette(cutree(fit, k = i), d)
    sil.width.val[counter] <- mean(clust[, 3])
  }

  # Combine the possible number of clusters and their silhouette widths into a data frame
  sil.res <- as.data.frame(cbind(sil.width.step, sil.width.val))

  # Determine the number of clusters to use (either based on the input or the silhouette widths)
  select.clst.num <- ifelse(is.null(part), sil.res$sil.width.step[sil.res$sil.width.val == max(sil.res$sil.width.val)], part)

  # Calculate the final silhouette data for the selected number of clusters
  final.sil.data <- cluster::silhouette(cutree(fit, k = select.clst.num), d)
  row.names(final.sil.data) <- group_labels

  # Get cluster assignments
  cluster_assignments <- cutree(fit, k = select.clst.num)

  # Generate colors for each individual boxplot, with colour defined by cluster
  cluster_colors <- rainbow(length(unique(cluster_assignments)), alpha = 1)[cluster_assignments]

  #create a vector of unit names ordered on the basis of how they appear in the cluster dendrogram
  labels.to.use <- group_labels[fit$order]

  # Create boxplot with colors according to the cluster assignment
  # on the x-axis, the boxplot will be given the ordered labels
  boxplot(groups[fit$order], col = cluster_colors[fit$order], names=labels.to.use, main="Boxplots colored by cluster membership", cex.main=0.90, cex.axis=0.70)

  # Plot the dendrogram
  graphics::plot(fit,
                 main = "Clusters Dendrogram",
                 sub = paste0("\nAgglomeration method: ", aggl.meth),
                 xlab = "",
                 cex = cex.dndr.lab,
                 cex.main = 0.90,
                 cex.sub = 0.75)

  # Add rectangles to the dendrogram to visualize the clusters
  solution <- rect.hclust(fit, k = select.clst.num, border = unique(cluster_colors[fit$order]))

  # Conditonally define the by-cluster colours to be used when
  # plotting the silhouette plot
  if(silh.col==TRUE){
    by.clust.col <- rainbow(length(unique(cluster_assignments)), alpha = 1)
  } else {
    by.clust.col <- NULL
  }

  # Plot the silhouette plot
  graphics::plot(final.sil.data,
                 col=by.clust.col,
                 cex.names = cex.sil.lab,
                 max.strlen = 30,
                 nmax.lab = num_groups + 1,
                 main = "Silhouette plot",
                 cex.main = 0.95)

  # Add a vertical line to show the average silhouette width
  abline(v = mean(final.sil.data[, 3]), lty = 2)

  # Plot the average silhouette width vs. the number of clusters
  graphics::plot(sil.res, xlab = "Number of clusters",
                 ylab = "Average silhouette width",
                 ylim = c(0, 1),
                 xaxt = "n",
                 type = "b",
                 main = "Average silhouette width vs. number of clusters",
                 sub = paste0("values on the y-axis represent the average silhouette width at each cluster solution"),
                 cex.main = 0.9,
                 cex.sub = 0.75)

  # Add an x-axis to the plot
  axis(1, at = 0:max.ncl, cex.axis = 0.7)

  # Add a point to show the selected number of clusters
  graphics::points(x = select.clst.num, y = sil.res[select.clst.num - 1, 2], pch = 20)

  # Add a label to the point with the silhouette width value
  graphics::text(x = select.clst.num, y = sil.res[select.clst.num - 1, 2], labels = round(sil.res[select.clst.num - 1, 2], 3), cex = 0.65, pos = 3, offset = 1.2, srt = 90)

  # Rename the columns of the sil.res dataframe to give
  # more meaningful labels before returning it
  colnames(sil.res)[1] <- "cluster solution"
  colnames(sil.res)[2] <- "average silhouette width"

  return(list("distance.matrix"=d,
              "units.by.cluster"=split(groups, cluster_assignments),
              "avr.silh.width.by.n.of.clusters"=sil.res,
              "partition.silh.data"=final.sil.data))
}
