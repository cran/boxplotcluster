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
#' @param x A dataframe or matrix where each column represents a data series to be clustered.
#' @param calc.type A string specifying the units to be clustered (either "columns" or "rows"; the former is the default).
#' @param aggl.meth A string specifying the agglomeration method to be used in hierarchical
#' clustering. Defaults to "ward.D2". For other methods see \code{\link[stats]{hclust}}.
#' @param part An optional integer specifying the desired number of clusters. If not provided,
#' the function will select the number of clusters that maximizes the average silhouette width.
#' @param cex.dndr.lab A numeric specifying the character expansion factor for the labels in the
#' dendrogram plot. Defaults to 0.75.
#' @param cex.sil.lab A numeric specifying the character expansion factor for the labels in the
#' silhouette plot. Defaults to 0.75.
#' @param oneplot TRUE (default) or FALSE if the user wants or does not want the plots to be visualized
#' in a single window.
#'
#' @details The function first calculates the pairwise distance between each row or column of the input
#' dataset using the Ichino-Yaguchi dissimilarity measure (equations 7 and 8 in Arroyo-Maté-Roque (2006)).
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
#' The distance matrix is then used to perform hierarchical clustering. Also, the function calculates the
#' silhouette width for different numbers of clusters and selects the number of clusters
#' that maximizes the average silhouette width (unless a specific number of clusters is provided by the user).\cr
#'
#' The silhouette method allows to measure how 'good' is the selected cluster solution. If the parameter \code{part}
#' is left empty (default), an optimal cluster solution is obtained. The optimal partition is selected via an iterative
#' procedure which identifies at which cluster solution the highest average silhouette width is achieved.
#' If a user-defined partition is needed, the user can input the desired number of clusters using the parameter \code{part}.
#' In either case, an additional plot is returned besides the cluster dendrogram and the silhouette plot; it displays a
#' scatterplot in which the cluster solution (x-axis) is plotted against the average silhouette width (y-axis).
#' A black dot represents the partition selected either by the iterative procedure or by the user.\cr
#'
#' In summary, the function generates a series of plots to visualize the results:\cr
#'
#' (a) boxplots colored by cluster membership,\cr
#' (b) a dendrogram (where clusters are indicated by rectangles whose color is consistent with the color assigned to the
#' boxplot in the previous plot),\cr
#' (c) a silhouette plot, and\cr
#' (d) a plot of the average silhouette width vs. the number of clusters.\cr
#'
#' The silhouette plot is obtained from the \code{silhouette()} function out from the \code{cluster} package.
#' For a detailed description of the silhouette plot, its rationale, and its interpretation, see Rousseeuw 1987.
#'
#' The function returns a list storing the following components \itemize{
##'  \item{distance.matrix: }{distance matrix reporting the distance values.}
##'  \item{dataset.w.cluster.assignment: }{a copy of the input dataset; if the rows are clustered, a new column is added
##'  which stored the cluster membership; if columns are clustered, a new row is added at the very end of the dataset
##'  to store the cluster membership.}
##'  \item{avr.silh.width.by.n.of.clusters: }{average silhouette width by number of clusters.}
##'  \item{partition.silh.data: }{silhouette data for the selected partition.}
##' }
#'
#' @references Rousseeuw, P J. (1987). Silhouettes: A graphical aid to the interpretation and validation of cluster analysis,
#' Journal of Computational and Applied Mathematics 20, 53-65.
#'
#' @references Ichino, M., & Yaguchi, H. (1994). Generalized Minkowski Metrics for Mixed
#' Feature-Type Data Analysis. IEEE Transactions on Systems, Man, and Cybernetics, 24(4), 698-708.
#'
#' @references Arroyo, J., Maté, C., & Roque, A. M-S. (2006). Hierarchical Clustering for Boxplot Variables.
#' In Studies in Classification, Data Analysis, and Knowledge Organization (pp. 59–66).
#' Springer Berlin Heidelberg.
#'
#'
#' @examples
#'
#'   # EXAMPLE 1
#'   # Create a toy dataset
#'   df <- data.frame(
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
#' # EXAMPLE 2
#' # Create a toy dataset representing archaeological stone flake length (cm) by raw material
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
#' boxplotcluster(df)
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
boxplotcluster <- function(x, calc.type = "columns", aggl.meth = "ward.D2", part = NULL, cex.dndr.lab = 0.75, cex.sil.lab = 0.75, oneplot=TRUE) {
  # Save current par settings
  oldpar <- par(no.readonly = TRUE)
  # Ensure settings are restored when function exits
  on.exit(par(oldpar))

  if(oneplot==TRUE){
    m <- rbind(c(1,2), c(3,4))
    layout(m)
  }

  # Transpose the matrix so rows become columns (and vice versa)
  if (calc.type == "rows") {
    x <- t(x)
  }

  # Function to calculate the distance between two columns of x
  calc_dist <- function(s1, s2) {
    m1 <- min(x[, s1])
    m2 <- min(x[, s2])
    q1 <- quantile(x[, s1], probs = 0.25)
    q2 <- quantile(x[, s2], probs = 0.25)
    Q1 <- quantile(x[, s1], probs = 0.75)
    Q2 <- quantile(x[, s2], probs = 0.75)
    M1 <- max(x[, s1])
    M2 <- max(x[, s2])
    Me1 <- median(x[, s1])
    Me2 <- median(x[, s2])

    (0.5 * (abs(m1 - m2) +
              2 * abs(q1 - q2) +
              2 * abs(Me1 - Me2) +
              2 * abs(Q1 - Q2) +
              abs(M1 - M2))) /4
  }

  # Calculate the distance matrix using the outer function and vectorized calc_dist function
  d <- outer(1:ncol(x), 1:ncol(x), Vectorize(calc_dist))

  # Convert the distance matrix to a proper labeled matrix
  d <- as.matrix(d, labels = TRUE)
  colnames(d) <- rownames(d) <- colnames(x)

  # Convert the distance matrix to a dist object
  d <- as.dist(d)

  # Perform hierarchical clustering
  fit <- hclust(d, method = aggl.meth)

  # Determine the maximum number of clusters
  max.ncl <- ncol(x) - 1

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
  row.names(final.sil.data) <- colnames(x)

  # Get cluster assignments
  cluster_assignments <- cutree(fit, k = select.clst.num)

  # Generate colors for each cluster
  cluster_colors <- rainbow(length(unique(cluster_assignments)))[cluster_assignments]

  # If the input dataframe (regardless of whether transposed or not) has no column labels,
  # the column numbers in the order appearing in the cluster dendrogram will be used as labels
  if(is.null(colnames(x)) == TRUE) {
    labels.to.use <- fit$order
  } else {
    # otherwise, the dataframe's column labels, in the order appearing in the cluster dendrogram, will be used
    # as labels
    labels.to.use <- colnames(x)[fit$order]
  }

  # Create boxplot with colors according to the cluster assignment
  # on the x-axis, the boxplot will be given the labels as defined in the previous if/else condition
  boxplot(x[,fit$order], col = cluster_colors[fit$order], names=labels.to.use, main="Boxplots colored by cluster membership", cex.main=0.90, cex.axis=0.70)

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

  # Plot the silhouette plot
  graphics::plot(final.sil.data,
                 cex.names = cex.sil.lab,
                 max.strlen = 30,
                 nmax.lab = ncol(x) + 1,
                 main = "Silhouette plot")

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

  # Create a copy of the input dataset and create either a new column or a new row to store the cluster membership
  # for each column/row
  if (calc.type == "rows") {
    x.copy <- as.data.frame(t(x))
    x.copy$cluster_assignments <- cluster_assignments
  } else {
    x.copy <- as.data.frame((x))
    x.copy[nrow(x.copy)+1,] <- round(as.numeric(cluster_assignments,0))
    rownames(x.copy[nrow(x.copy),]) <- "cluster assignment"
  }

  # Rename the columns of the sil.res dataframe to give
  # more meaningful labels before returning it
  colnames(sil.res)[1] <- "cluster solution"
  colnames(sil.res)[2] <- "average silhouette width"

  return(list("distance.matrix"=d,
              "dataset.w.cluster.assignment"=x.copy,
              "avr.silh.width.by.n.of.clusters"=sil.res,
              "partition.silh.data"=final.sil.data))
}
