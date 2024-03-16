---
title: "`spillR`: Spillover Compensation in Mass Cytometry Data"
author: "Marco Guazzini$^{1}$, Alexander G. Reisach$^{2}$, Sebastian Weichwald$^{3}$, and Christof Seiler$^{1,4,5}$"
date: "$^1$Department of Advanced Computing Sciences, Maastricht University, The Netherlands \\\n $^2$Université Paris Cité, CNRS, MAP5, F-75006 Paris, France \\\n $^3$Department of Mathematical Sciences, University of Copenhagen, Denmark \\\n $^4$Mathematics Centre Maastricht, Maastricht University, The Netherlands \\\n $^5$Center of Experimental Rheumatology, Department of Rheumatology, \\\n University Hospital Zurich, University of Zurich, Switzerland \\\n \\\n March 16, 2024"
output:
  bookdown::pdf_document2: 
    toc: false
    extra_dependencies: ["euflag"]
bibliography: spillr_paper.bib
csl: style.csl
link-citations: true
abstract: |
  Channel interference in mass cytometry can cause spillover and may result in miscounting of protein markers. @catalyst introduce an experimental and computational procedure to estimate and compensate for spillover implemented in their R package `CATALYST`. They assume spillover can be described by a spillover matrix that encodes the ratio between unstained and stained channels. They estimate the spillover matrix from experiments with beads. We propose to skip the matrix estimation step and work directly with the full bead distributions. We develop a nonparametric finite mixture model and use the mixture components to estimate the probability of spillover. Spillover correction is often a pre-processing step followed by downstream analyses, and choosing a flexible model reduces the chance of introducing biases that can propagate downstream. We implement our method in an R package `spillR` using expectation-maximization to fit the mixture model. We test our method on synthetic and real data from `CATALYST`. We find that our method compensates low counts accurately, does not introduce negative counts, avoids overcompensating high counts, and preserves correlations between markers that may be biologically meaningful.
---





# Introduction

Mass cytometry makes it possible to count a large number of proteins simultaneously on individual cells [@bandura2009mass; @bendall2011single]. Although mass cytometry has less spillover---measurements from one channel overlap less with those of another---than flow cytometry [@sp-c; @novo2013generalized], spillover is still a problem and affects downstream analyses such as differential testing [@diffcyt; @seiler2021cytoglmm] or dimensionality reduction [@scater]. Reducing spillover by careful design of experiment is possible [@takahashi2017mass], but a purely experimental approach may be neither efficient nor sufficient [@lun2017influence].

@catalyst propose a method for addressing spillover by conducting an experiment on beads. This experiment measures spillover by staining each bead with a single type of antibody. The slope of the regression line between target antibodies and non-target antibodies represents the spillover proportion between channels. @miao2021ab attempt to solve spillover by fitting a mixture model. Our contribution combines the solutions of @catalyst and @miao2021ab. We still require a bead experiment, as in @catalyst, but estimate spillover leveraging a statistical model, as in @miao2021ab. Both previous versions rely on an estimate for the spillover matrix. The spillover matrix encodes the pairwise spillover proportion between channels. We avoid estimating a spillover matrix and instead model spillover by fitting a mixture model to the observed counts. Our main new assumption is that the spillover distribution---not just the spillover proportion---from the bead experiment carries over to the biological experiment. In other words, we transfer the spillover distribution to the real experiment instead of just the spillover proportion encoded in the spillover matrix.

We present our mixture model and link it to calculating spillover probabilities for specific count values in Section \@ref(methods). Our estimation procedure is based on an EM algorithm and logistic regression, and implemented in our new R package `spillR`^[\textcolor{red}{https://bioconductor.org/packages/spillR}]. We conduct experiments on simulated, \textcolor{red}{semi-simulated}, and real data obtained from the `CATALYST` R package [@catalyst] in Section \@ref(results), and discuss our experiments and relate our findings to `CATALYST` in Section \@ref(discussion).

# Methods

\textcolor{red}{
In this section we first illustrate our method \texttt{spillR} (as well as a simple baseline version \texttt{spillR (naive)} by an example, and then describe the algorithm as well as its underlying assumptions. Regarding the terminology, note that in mass cytometry, counts are often referred to as dual counts or signal intensity. We call them counts to emphasize that we rely on the fact that they are non-negative integers as opposed to possibly real valued intensities. 
}

## Example



\begin{figure}

{\centering \includegraphics[width=0.75\linewidth]{spillr_paper_files/figure-latex/method-example-1} 

}

\caption{\textcolor{red}{Panel A shows a density plot of target and spillover markers, Panel B shows spillover probability for Yb173Di estimated by \texttt{spillR}, and Panel C compares spillover compensation by our methods and \texttt{CATALYST}. Counts are arcsinh transformed with cofactor of five (\protect\hyperlink{ref-bendall2011single}{Bendall \emph{et al.}, 2011}), zero counts are not shown; they are 31 for no compensation, 1603 for \texttt{spillR}, 2162 for \texttt{CATALYST}, and 2 nbins for \texttt{spillR (naive)}. As shown in Panel C, our baseline method \texttt{spillR (naive)} performs similarly to \texttt{CATALYST} and removes the first peak of the uncorrected data (red) between about 2 and 4. By contrast, \texttt{spillR} is sensitive to the difference in shape between the peaks in the bead data (Panel A) the first peak in the real data (Panec C red), and only removes the part that corresponds to a peak in the bead experiment.}}(\#fig:method-example)
\end{figure}

Our procedure is illustrated in Figure \@ref(fig:method-example), using a dataset from the `CATALYST` package as an example. There are four markers, HLA-DR (Yb171Di), HLA-ABC (Yb172Di), CD8 (Yb174Di), and CD45 (Yb176Di), that spill over into the target marker, CD3 (Yb173Di). The markers have two names: the first name is the protein name and the second name in brackets is the conjugated metal. There are bead experiments for each of the spillover markers.

Panel A depicts the marker distributions from the beads experiment. We see that for this marker the bead experiments are high-quality as the target marker Yb173Di is concentrated around six, similarly to the experiment with real cells. This suggests that the spillover marker values can be transferred to the real experiments. Marker Yb172Di shows large spillover into Yb173Di, and suggests that the left tail of the first mode of the distribution may be attributed to that marker. The other spillover markers have low counts, making it justifiable to set some or all the low counts to zero.

\textcolor{red}{Panel B shows a curve representing our spillover probability estimates. We can see that the probability of spillover is highest at points that correspond to a high density of spillover markers in Panel A. If the spillover probability is close to one, our correction step assigns most cells to spillover.} Counts above four stem from spillover with probability zero (and from the actual target with probability one), which means that our procedure keeps them at their raw uncorrected value.



Panel C displays the distribution of our target marker, CD3 (Yb173Di) before and after spillover correction. \textcolor{red}{We observe few real counts (red) below a value of $2$, so although all methods perform strong compensation in this range, there is little visible change. From $2$ onward there is a clear distinction between the methods. `CATALYST`, like our baseline `spillR (naive)`, compensates nearly all counts below the second peak of the raw counts (red) as spillover. By contrast, `spillR` compensates only where the relative frequency of the raw counts (red) match the density of spillovermarkers in the bead experiment shown in panel A, and does not remove a part of the first peak as a result. While `CATALYST` shifts the distribution of large counts (around 6) slightly to the left, our methods leave them unaffected as the bead experiment shows no spillover in this region.}
<!-- Zero counts (or equivalently `NA` counts for `spillR`) and mean counts are shown in Table \@ref(tab:method-example-summary). 
As we do not know the true mean in this case, the comparison of means here serves merely as an illustration of the main differences between \texttt{spillR} and \texttt{CATALYST} and not to identify the correct method.
-->
\textcolor{red}{
Our baseline method \texttt{spillR (naive)} is similar to \texttt{CATALYST} in low and medium range, but keeps higher counts unchanged.
}

## Definition of Spillover Probability and Assumptions

We observe a count $Y_i$ of a target marker in cell $i$. We model the observed $Y_i$ as a finite mixture [@mclachlan2019finite] of unobserved true marker counts $Y_i \mid Z_i = 1$ and spillover marker counts $Y_i \mid Z_i = 2, \dots, Y_i \mid Z_i = K$ with mixing probabilities $\pi_{k} = P(Z_i = k)$ for $k = 1, \dots, K$, 
$$
P(Y_i = y) = \sum_{k = 1}^K \pi_k \, P(Y_i = y \mid Z_i = k).
$$
The first mixing probability is the proportion of true signal in the observed counts. The other $K-1$ mixing probabilities are the proportions of spillover. The total sum of mixing probabilities equals one, $\sum_k \pi_k = 1$. The total number of markers in mass cytometry panels is between 30 and 40 [@bendall2011single], but only a small subset of three to four markers spill over into the target marker [@catalyst]. So, typically $K = 1+3$ or $K = 1+4$.

Experimentally, we only measure a sample from the distribution of $Y_i$. The probabilities $\pi_k$ and true distributions $P(Y_i = y \mid Z_i = k)$ are unobserved, and we need to estimate them from data. 
In many applications, the mixture components are \textcolor{red}{modeled to be} in a parametric family, for example, the negative binomial distribution.
As spillover correction is a pre-processing step followed by downstream analyses, choosing the wrong model can introduce biases in the next analysis step. To mitigate such biases, we propose to fit nonparametric mixture components. We make two assumptions that render the components and mixture probabilities identifiable:

* <div id='assumption1'>(A1) Spillover distributions are the same in bead and real experiments.</div> 
The distribution of $Y_i \mid Z_i = k$ for all $k > 1$ is the same in beads and real cells. This assumption allows us to learn the spillover distributions of $Y_i \mid Z_i = k$ for all $k > 1$ from experiments with beads, and transfer them to the experiment with real cells. This assumption relies on high-quality single-stained bead experiments that measure spillover in the same range as the target biological experiment. In other words, a high-quality bead experiment for our method works best if the distribution of bead cells is similar to the distribution of real cells.

* <div id='assumption2'>(A2) For each cell $i$, the observed count $Y_i$ can only be due to one distribution.</div> 
This assumption is already implied by the statement of the mixture model. It allows us to calculate the spillover probability for a given count $Y_i = y$ from the posterior probability that it arises through spillover from markers $k > 1$,
$$
P(\text{spillover} \mid Y_i = y) = P(Z_i > 1 \mid Y_i = y) = 
1 - P(Z_i = 1 \mid Y_i = y) = 
1 - \frac{\pi_1 \, P(Y_i = y \mid Z_i = 1)}{P(Y_i = y)}.
$$
To parse this calculation, recall that in mixture models the $\pi_1$ is the prior probability, $P(Y_i = y \mid Z_i = 1)$ is the conditional probability given the mixture component, and the denominator $P(Y_i = y)$ is the marginal distribution. Applying Bayes rule leads to the posterior probability.

## Estimation of Spillover Probability

We propose a two step procedure for estimating the spillover probability. In step 1, we estimate mixture components and mixture probabilities. We refine these estimates using the EM algorithm [@dempster1977maximum]. In step 2, we use these probability estimates to assign counts to spillover or signal.

We denote the $n \times K$ count matrix as $\mathbf{Y} = (y_{ik})$ with real cells in the first column and beads in columns two and higher. To simplify mathematical notation but without loss of generality, we assume that the number of 
\textcolor{red}{events}
from real and bead experiments have the same $n$. In practice, the number of 
\textcolor{red}{events}
from bead experiments is much smaller than from real experiments. The $k$th column of $\mathbf{Y}$ contains marker counts for the $k$th spillover marker, which represents the empirical spillover distribution of marker $k$ into the target marker, that is, the marker in the first column of $\mathbf{Y}$.

### EM Algorithm

* Initialization: For the mixture probability vector, we assign probability $0.9$ to the the target marker and divide the probability $0.1$ among the spillover markers, 
$$
\hat{\pi}_{1} = 0.9 \text{ and } \hat{\pi}_i = 0.1/(K-1) \text{ for all } i > 1.
$$ 
The procedure is not sensitive to the choice of the initial mixture probability vector and other initializations are possible but may be slower to converge. Then, we initialize the $k$th mixture component using \textcolor{red}{its probability mass function (PMF) after smoothing and normalizing, $\widehat{P}(Y_i = y \mid Z_i = k)$. We smooth the PMF using kernel density estimation implemented in R function \texttt{density} with the default option for selecting the bandwith of a Gaussian kernel.}

* E-step: We evaluate the posterior probability of a count $y$ belonging to component $k$ (that is, originating from marker $k$),
$$
\widehat{P}\left(Z_i = k \mid Y_i = y \right) = 
\frac
{ \hat{\pi}_k \, \widehat{P}(Y_i = y \mid Z_i = k) }
{ \sum_{k' = 1}^K \hat{\pi}_{k'} \, \widehat{P}(Y_i = y \mid Z_i = k') }.
$$

* M-step: We estimate the new mixture probability vector from posterior probabilities, 
$$
\hat{\pi}_k = 
\frac{1}{n} \sum_{i = 1}^n \widehat{P}\left(Z_i = k \mid Y_i = y \right),
$$
and estimate the new target marker distribution \textcolor{red}{by smoothing and normalizing. Here, we use the R function \texttt{density} again}, but weight each observation according to their posterior probabilities, $\widehat{P} \left(Z_i = 1 \mid Y_i = y \right)$. We only update the target marker, $\widehat{P}(Y_i = y \mid Z_i = 1)$, and keep the other bead distributions, $\widehat{P}(Y_i = y \mid Z_i = k)$ for all $k > 1$, fixed at their initial value.

To refine our estimates, we iterate over the E and M-steps until estimates stabilize. We stop iterating when $\hat{\pi}_1$ changes less than $10^{-5}$ from the previous iteration. The final output is the spillover probability curve with estimates at discrete points in the support of $Y_i$, 
$$
\widehat{P}(\text{spillover} \mid Y_i = y) = 1 - \widehat{P}(Z_i = 1 \mid Y_i = y).
$$

We rely on assumption [(A1)](#assumption1) to justify updating only the distribution of the target marker. We rely on assumption [(A2)](#assumption2) to justify calculating the spillover probability from the mixture model. We refer to Appendix \@ref(em-algorithm-example) for a step-by-step example of our EM algorithm. 

### Spillover Decision

To perform spillover compensation, we draw from a Bernoulli distribution with the spillover probability as parameter to decide whether or not to assign a given count to spillover. \textcolor{red}{We set counts designated to spillover to a user-specified value. We recommend a value of zero to maintain the overall cellular composition of the sample or a value such as \texttt{NA} or $-1$ to mark them explicitly as spillover for downstream processing.}

<!-- 

## Identifiability

[TODO] Need to think if our model assumptions guarantee identifiability and if our procedure is guaranteed to converge. I think if it is identifiable, then convergence follows from the property of the EM algorithm. Some related work [@hall2003nonparametric; @aragam2020identifiability] that might help. 

-->

## Baseline Method `spillR (naive)`

\textcolor{red}{
We compare our mixture method to a naive baseline method that considers only the bead distributions. Similarly to our standard \texttt{spillR} method, we estimate the bead PMF of each bead $k$ with the kernel density estimator \texttt{density}, $\widehat{P}(Y_i | Z_i = k)$. Then, for all count values $y$ in the range of the bead counts, we separately normalize the PMF at each value $Y_i = y$ and calculate the spillover probability as, 
$$
P(Z_i = k \mid Y_i = y) = 1 - \frac{\widehat{P}(Y_i | Z_i = k)}{\sum_{k'}^K \widehat{P}(Y_i | Z_i = k)}.
$$
The spillover decision is the same as in our standard \texttt{spillR} method. This is a computationally efficient and simple baseline that, in effect, assigns counts to the marker with the highest density at the count value in the corresponding bead experiment.
}

# Results

We first evaluate our new method `spillR` on simulated datasets. We probe our method to experimentally find its shortcomings. Then, we compare `spillR` to the non-negative least squares method implemented in the R package `CATALYST` on real and \textcolor{red}{semi-simulated} data from the same package. All experiments and plots can be reproduced by compiling the R markdown file `spillR_paper.Rmd`^[https://github.com/ChristofSeiler/spillR_paper].

## Simulated Data





\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{spillr_paper_files/figure-latex/simulated-experiments-plot-1} 

}

\caption{Three experiments testing our assumptions and sensitivity to bimodal bead distribution. For each experiment the top row are mean values over the entire range of the experimental setups. \textcolor{red}{The mean values for \texttt{spillR} are computed with \texttt{NA} imputation for spillover counts, so the mean is identical to the true mean without spillover if all spillover counts are correctly identified as such.} The bottom row are density plots for three parameter settings to illustrate the generated distributions. $Y$ is the distribution with spillover. $Y \mid Z = 1$ is the distribution without spillover. $Y \mid Z = 2$ is the spillover. mean($Y$) is the average of the distribution with spillover. mean($Y \mid Z = 1$) is the average count without spillover. \texttt{spillR} mean($Y$) is the average count after correcting $Y$.}(\#fig:simulated-experiments-plot)
\end{figure}

We choose three different experiments to test `spillR` against different bead and real cell distributions. We explore a wide range of possible parameter settings. Figure \@ref(fig:simulated-experiments-plot) has three panels, each representing one experimental setup. The first two panels test our assumptions [(A1)](#assumption1) and [(A2)](#assumption2). The third panel tests sensitivity of `spillR` to bimodal bead distributions. For all three experiments, we model counts using a Poisson distribution with parameter $\lambda$. We simulate 10,000 real cells with $\lambda = 200$, and 1,000 beads with $\lambda = 70$, and a spillover probability of $0.5$. The bead data are an independent copy of the true spillover. The other parameters and statistical dependencies are specific to each experiment. The details of the generative models are given in Appendix \@ref(generative-models). We repeat each simulation 20 times and report averages over the 20 replications. 

Each panel of Figure \@ref(fig:simulated-experiments-plot) has two rows of plots. The plot in the first row represents the summary of the means for each experimental setup as a function of their respective parameter $\tau$. This parameter has a different meaning in each setup. To visualize the different experiments, we summarize the full distributions with the true simulated signal mean (black), the uncorrected mean (orange), and the `spillR` corrected mean (green). Plots on the second row illustrate the simulated data distributions for three selected parameters $\tau$ picked from the experimental setup. The yellow density curve shows the observed counts $Y$. The black density curve shows the distribution of target cell counts. The blue density curve shows the distribution of spillover counts. The goal of the experiment is to estimate the mean of the black density as accurately as possible from the yellow density curve, which represents the data $Y$ we would observe in practice. We simulate this data using the models in Appendix \@ref(generative-models).

<!--
### Bead Shift (A1)
-->

In the first experiment (panel A), we shift the spillover in the beads experiment away from the true spillover to probe [(A1)](#assumption1). We test a range of bead shifts from no shift at $\tau = 0$ to \textcolor{red}{$\tau = -10$. At $\tau = -10$,} the measured spillover (the first mode of the yellow density) is shifted away from the actual spillover (the blue density), causing both the observed and compensated mean to be lower than the true mean. This may be the case in a low-quality bead experiment. As $\tau$ gets closer to zero, the first mode of the yellow density moves towards the blue density (as may be the case in a higher quality bead experiment), and the compensated signal moves closer to the true mean.
<!-- Our compensation also increases means, as it should, in contrast to e.g. @catalyst. -->

<!--
### Model Misspecification (A2)
-->

In the second experiment (panel B), we mix target and spillover to explore the robustness of our method with respect to our second assumption [(A2)](#assumption2). One way to think about this is that the mixture is a form of model misspecification. Our mixture model is undercomplete, which means that there are more true mixture components than we observe in the beads experiment. If $\tau = 0$, then assumption [(A2)](#assumption2) is correct, but for $\tau = 0.5$ the assumption [(A2)](#assumption2) is maximally violated. The true mean decrease with increasing $\tau$. 
\textcolor{red}{\texttt{spillR} compensates well as long as $\tau$ is close to zero, but deviates from the true mean with increasing $\tau$. As the spillover distribution becomes more similar to the target marker distribution, the mean of the \texttt{spillR} flips to the mean of the observed data at $\tau \approx 0.25$, until}
<!-- Our compensation is closer to the true mean across the tested range. -->
at $\tau = 0.5$ all three distributions and their means are the same.

<!--
### Bimodal Spillover
-->

In the third experiment (panel C), we model spillover with a bimodal distribution. Here $\tau$ is the mixing probability of the two modes. The locations of the two spillover modes are fixed. If $\tau = 0$ or $\tau = 1$, then spillover is unimodal. If $\tau = 0.5$, the first mode of the bimodal beads distribution is left to the signal mode, and the second mode is to the right. The corrected mean is closer to the true mean than the uncorrected mean across the test range.

## Real Data



\begin{figure}

{\centering \includegraphics[width=0.95\linewidth]{spillr_paper_files/figure-latex/spillr-vignette-1} 

}

\caption{Comparison of compensation methods and uncorrected counts on real data. Counts are arcsinh transformed with cofactor of five (\protect\hyperlink{ref-bendall2011single}{Bendall \emph{et al.}, 2011}).}(\#fig:spillr-vignette)
\end{figure}


We compare our methods to `CATALYST` on one of the example datasets in the `CATALYST` package. The dataset consists of an experiment with real cells and corresponding single-stained bead experiments. The experiment on real cells has 5,000 peripheral blood mononuclear cells from healthy donors measured on 39 channels. The experiment on beads has 10,000 cells measured on 36 channels with the number of beads per metal label ranging from 112 to 241.

In Figure \@ref(fig:spillr-vignette), we show the comparison of our methods to `CATALYST` on the same markers as their original paper [@catalyst] in their Figure 3B. In the original experiment, they conjugate the three proteins CD3, CD8, and HLA-DR with two different metal labels. 
\textcolor{red}{
For example, they conjugate CD8 with Yb174Di (Yb is the metal and the number indicates the number of nucleons of the isotope) and La139Di. 
As in their plot, our columns correspond to the different metal labels.
Following their set-up, we show the three target proteins on the vertical axis.
On the horizontal axis we show the spillover markers; we show CD3 in row one, and HLA-ABC in row two and three.
}
We visualize the joint distributions using two-dimensional histograms.

In all six panels (A--F), we observe that `spillR` compensates most strongly in the low counts. In panel C, CD3 (Yb173Di) against HLA-ABC (Yb172Di), `CATALYST` can be seen to compensate strongly in the middle range. It removes the spherical pattern that shows correlation between the two markers. `spillR` preserves this correlation structure and only masks out the lower counts of CD3 (Yb173Di). This highlights a key difference between `spillR` and `CATALYST`: `spillR` \textcolor{red}{identifies counts that may arise from spillover and replaces them with a user-specified value (e.g. 0, \texttt{NA}, or -1), whereas `CATALYST` shrinks counts across the entire range to compensate for spillover.}
\textcolor{red}{
\texttt{spillR (naive)} compensates most aggressively, assigning most counts outside the target marker bead distribution to spillover.}


\textcolor{red}{
In Panels D and F we can see that \texttt{CATALYST} performs no visible compensation. For the same panels our \texttt{spillR} methods compensate strongly. Our diagnostic plots using function \texttt{plotDiagnostics} in \texttt{spillR} give some indication on why our method compensates that strongly: The original experiment is designed such that we expect no spillover in the channels on the vertical axis, but our diagnostic plots show that the bead distribution of at least one spillover marker---as identified by the spillover matrix in \texttt{CATALYST}---overlaps with the first mode in the distribution of real cells. Thus, the strong compensation by \texttt{spillR} is consistent with our assumption (A1). The results of \texttt{spillR (naive)} fit with this explanation as it completely removes all counts in that area, which occurs when there is only spillover and no signal from the target marker.
} 
<!-- [(A1)](#assumption1) -->

The color code of the two-dimensional histograms indicates the absolute number of cells that fall into one hexagon bin. The uncorrected and `spillR` corrected histograms can contain different absolute numbers of cells,
\textcolor{red}{
even for identical distributions due to a rounding step in \texttt{spillR} that converts raw counts to integers. Raw mass cytometry data may not be true count data because the proprietary post-processing of the manufacturer often performs a randomization step when exporting the data.
}
The uncorrected counts do not undergo this pre-processing step. `CATALYST` does not perform this pre-processing step. This also explains the different patterns in panel B. `spillR` has horizontal stripes that correspond to non-integer values not in the support of the distribution for `spillR`.
\textcolor{red}{
We leave the decision to apply re-randomization of the count data for downstream analysis up to the user. Our rational is that the user should see the differences in this pre-processing step and how it propagates to the results.
}



\textcolor{red}{
The average computation time with 100 replications on an Apple M1 with 8 cores and 16 GB of RAM is $10.6$ seconds for \texttt{spillR}, $0.43$ seconds for \texttt{CATALYST}, and $0.45$ seconds for \texttt{spillR (naive)}. The computational costs scale linearly in the number of cells and number of spillover markers. This allows for processing large-scale datasets.
}

## Semi-Simulated Data



\begin{figure}

{\centering \includegraphics[width=1\linewidth]{spillr_paper_files/figure-latex/semi-simulated-plot-1} 

}

\caption{Comparison of compensation methods and uncorrected counts on semi-synthetic data 	extcolor{red}{(	exttt{spillR} and 	exttt{spillR (naive)} are set to impute spillover values with $0$)}. The vertical dashed line helps to interpret the spillover correction. It indicates the 	extcolor{red}{original mode of the bead distribution of Yb172Di at 2.725, before it overwriting it with the first peak of the real observations of Yb172Di}. Counts are arcsinh transformed with cofactor of five (\protect\hyperlink{ref-bendall2011single}{Bendall \emph{et al.}, 2011}).}(\#fig:semi-simulated-plot)
\end{figure}



\textcolor{red}{
We compare \texttt{spillR} and \texttt{CATALYST} on semi-simulated data. The goal is to elucidate differences between \texttt{spillR} and \texttt{CATALYST}, and to evaluate the performance of \texttt{spillR} if more than one marker spills into the target marker. We create semi-simulated datasets by overwriting the beads distribution for the target marker CD3 (Yb173Di). We consider the first mode of the count distribution of CD3 (Yb173Di) observed in real cells, that is, the counts from $1.44$ to $4.79$ on the transformed scale. We overwrite the beads distribution and shift counts in this range (mostly corresponding to Yb172Di) by three different values: no shift is $0$, subtracting $0.47$ on the transformed scale, and subtracting $0.94$ on the transformed scale. We further subsample without replacement from this new bead distribution to keep the same number of beads as in the original dataset. Figure \ref{fig:semi-simulated-plot} shows the three different beads experiment datasets in row A and the resulting compensations in row B.
}

\textcolor{red}{
In the first column of Figure \ref{fig:semi-simulated-plot}, the bead distributions are equal to the original dataset from Figure \ref{fig:method-example} except Yb172Di is now perfectly aligned with the first mode of the distribution of real cells (red curve in row B). In the second and third column, we shift the bead distribution of Yb172Di by $0.47$ and $0.94$. All three methods correctly compensate the spillover mode when no shift is present (first column). \texttt{CATALYST} and \texttt{spillR (naive)} compensates more aggressively in the medium shift cases (second column), while \texttt{spillR} is more moderate and compensates only the left hand tail of the spillover mode. For a shift of $0.94$ (third column), the three methods differ: \texttt{CATALYST} shrinks counts towards zero, shifting the entire spillover towards zero, \texttt{spillR} compensates lightly on the left hand tail, and \texttt{spillR (naive)} compensates aggressively leaving only a small right hand tail. This experiment illustrates how \texttt{spillR} compensates most strongly for counts that can be attributed to spillover following the distribution observed in the beads experiment.
}

# Discussion

\textcolor{red}{The sensitivity analysis in Section 3.1 illustrates the performance of \texttt{spillR} when our assumptions are violated.}
The experiment for [(A1)](#assumption1) shows that the mean count after `spillR` correction is closer to the true mean over a wide range of bead shifts. This indicates that our method can perform well even if the bead experiments are imperfect. If the difference between distributions of beads and real cells is large, then one option is to rerun the bead experiments to reduce this gap. The experiment for [(A2)](#assumption2) shows that our method is also robust to model misspecification. Additionally, misspecification can be addressed by adding all channels if necessary. The increase in computational cost when adding channels is relatively minor as our method scales linearly in the number of spillover markers. The experiment on bimodal bead distributions shows that the mean count after correction is still closer to the true mean even with bimodal bead distributions, even if the spillover is larger than the true signal.

In our comparison with `CATALYST` on \textcolor{red}{real data (Section 3.2) and semi-simulated data (Section 3.3)}, we observe the effect of the two different correction strategies. `CATALYST` essentially shrinks counts towards zero by minimizing a non-negative least squares objective. It assumes that spillover is linear up to counts of 5,000. The applied shrinkage is the same for low counts (e.g., below 10) and high counts (e.g., more than 100). By contrast, `spillR` does not require linearity of the spillover, but assumes that the distribution on the beads experiment carries over to the real cells experiment. In other words, the optimal beads experiment has the same peaks as the real cells experiment.

If counts are in the spillover range (which mostly applies to low counts), they are corrected strongly and set to \textcolor{red}{a user-specified imputation value}. If counts are not in the spillover range, they are left unchanged.
\textcolor{red}{Among the unchanged counts,} correlations between markers are preserved. The marker correlation between HLA-ABC (Yb172Di) and CD3 (Yb173Di) shown in the first column of Figure \@ref(fig:spillr-vignette) illustrates this point. `CATALYST` removes the positively correlated count concentration, whereas `spillR` keeps it. Compensation methods should try to remove spillover while keeping potentially biologically meaningful signal for unbiased downstream analyses. In this example, further experiments on the correlation structure between these markers would be necessary to resolve the discrepancy between the two methods. This is an important point as discovering correlations between markers can lead to the discovery of new clusters or signaling networks. 

\textcolor{red}{
Our baseline method, \texttt{spillR (naive)}, illustrates the behavior of more aggressive compensation by considering only the bead distribution. If our baseline method compensates aggressively in a certain range, this is because there are only non-target beads observed in that range. This approach highlights the allure and pitfalls of overcorrecting. While in Figure 1 it may seem that \texttt{spillR (naive)} removes all spillover by removing the first mode (Just like \texttt{CATALYST}), a closer inspection reveals that discrepancies between the bead spillover distribution and the first mode of the real cell distribution are not taken into account by either method, but do reflect in the compensation of \texttt{spillR}. A similar pattern can be seen in panels A, B, and C of Figure 3. This behavior reflects our assumption (A1) and highlights the role of bead experiments in the compensation process performed by \texttt{spillR}.
}

An additional advantage of our method is the diagnostic plot of the spillover probability curve. We can judge if the curve makes sense by comparing it to the observed count and bead distributions. Methods based on non-negative least squares \textcolor{red}{such as \texttt{CATALYST}} are harder to diagnose as they minimize a cost function with no clear biological interpretation. In our view, one of the biggest strengths of our current method is that it does not assume a specific parametric model for count data. We believe that this is crucial because spillover compensation precedes many downstream analysis steps, and avoiding the introduction of bias is thus our top priority.

\textcolor{red}{
For future work, we plan to apply our methodology to imaging mass cytometry (Angeloet al., 2014; Giesenet al., 2014; Bodenmiller, 2016). The basic method can be applied as is, but it will be beneficial to incorporate a spatial regularization term that enforces neighboring spillover to be similar to one another.
}

# Acknowledgments {-}

We thank EuroBioC2022 for awarding Marco Guazzini a travel award to present a preliminary version of `spillR` in Heidelberg. We thank Antoine Chambaz for his feedback on an earlier draft that substantially improved the paper. Alexander G. Reisach received funding from the European Union's Horizon 2020 research and innovation program under the Marie Sk\l{}odowska-Curie grant agreement No 945 nbins \euflag.

# References {-}

<div id="refs"></div>

\newpage

# (APPENDIX) Supplementary Material {-}

# EM Algorithm Example

Here we illustrate the procedure using a numerical example that includes one target and one spillover marker. We have one data matrix $\mathbf{Y}$ that contains real cell counts recorded for marker 1 (column 1) and the bead counts for marker 1 when the true marker was marker 2 (column 2). In practice, $\mathbf{Y}$ is usually a matrix with more than two columns representing multiple spillover markers. The index $i$ is a specific cell in beads and real cells experiment, respectively. Let's assume the following counts,
$$
\mathbf{Y} = (y_{ij}) = 
\begin{bmatrix}
3 & 2 \\ 
5 & 3 \\ 
17 & 2 \\ 
3   \\ 
17 \\ 
2 
\end{bmatrix}.
$$


```r
target    <- c(3, 5, 17,  3,  17, 2)
spillover <- c(2, 3,  2, NA, NA, NA)
Y = dplyr::bind_cols(target = target, spillover = spillover)
Y
```

```
## # A tibble: 6 x 2
##   target spillover
##    <dbl>     <dbl>
## 1      3         2
## 2      5         3
## 3     17         2
## 4      3        NA
## 5     17        NA
## 6      2        NA
```

* Initialization: We initialize our EM algorithm by estimating the conditional probability of observing $y$ given that it belongs to the target marker, and another conditional probability given that it belongs to the spillover marker.


```r
y_min <- min(Y$target)
y_max <- max(Y$target)
y_support <- y_min:y_max
fit1 <- density(Y$target, from = y_min, to = y_max)
fit2 <- density(Y$spillover, from = y_min, to = y_max, na.rm = TRUE)
f1 <- approxfun(fit1$x, fit1$y)
f2 <- approxfun(fit2$x, fit2$y)
P_Y1 <- f1(y_support)
P_Y1 <- P_Y1 / sum(P_Y1)
P_Y2 <- f2(y_support)
P_Y2 <- P_Y2 / sum(P_Y2)
P_YZ <- dplyr::bind_cols(P_Y1 = P_Y1, P_Y2 = P_Y2)
```

We initialize the mixture probabilities with the discrete uniform.


```r
pi <- c(0.9, 0.1)
```

Now, we update these initial values using the E and M-steps.

* E-step: Calculate the posterior probability for the true marker, and the spillover marker.


```r
P_ZY <- dplyr::mutate(P_YZ, 
                      P_Y1 = pi[1] * P_Y1, 
                      P_Y2 = pi[2] * P_Y2)
P_ZY <- P_ZY / rowSums(P_ZY)
P_ZY <- dplyr::bind_cols(target = y_support, P_ZY)
```

* M-step: Update the mixing probability vector,


```r
n <- nrow(Y)
YP <- dplyr::left_join(Y, P_ZY, by = "target")
YP
```

```
## # A tibble: 6 x 4
##   target spillover  P_Y1     P_Y2
##    <dbl>     <dbl> <dbl>    <dbl>
## 1      3         2 0.717 2.83e- 1
## 2      5         3 1.00  5.70e-13
## 3     17         2 1     9.29e-17
## 4      3        NA 0.717 2.83e- 1
## 5     17        NA 1     9.29e-17
## 6      2        NA 0.550 4.50e- 1
```

```r
pi <- c(sum(YP$P_Y1) / n, sum(YP$P_Y2) / n)
pi
```

```
## [1] 0.8307516 0.1692484
```

and re-estimate the distribution for the target marker using the posterior probabilities as weights, keep the non-target marker at its initial value,


```r
fit1 <- density(Y$target, from = y_min, to = y_max, weights = YP$P_Y1)
f1 <- approxfun(fit1$x, fit1$y)
P_Y1 <- f1(y_support)
P_Y1 <- P_Y1 / sum(P_Y1)
P_YZ <- bind_cols(P_Y1 = P_Y1, P_Y2 = P_Y2)
```

and calculate the spillover probability estimate,


```r
P_ZY <- dplyr::mutate(P_YZ, 
                      P_Y1 = pi[1] * P_Y1, 
                      P_Y2 = pi[2] * P_Y2)
P_ZY <- P_ZY / rowSums(P_ZY)
P_ZY <- dplyr::bind_cols(target = y_support, P_ZY)
P_ZY |>
  dplyr::mutate(p_spillover = round(1 - P_Y1, digits = 3)) |>
  dplyr::select(target, p_spillover) |>
  dplyr::filter(target %in% unique(Y$target))
```

```
##   target p_spillover
## 1      2       0.631
## 2      3       0.449
## 3      5       0.000
## 4     17       0.000
```

This is the result after one iteration.

# Generative Models

## Bead Shift {-}

Generative model for real cells $Y$ of this experiment:
$$
\begin{aligned}
I              & \sim \text{Bernoulli}(0.1)                            & \qquad \text{(spillover indicator)} \\
Z              & = I + 1                                               & \qquad \text{(channel number)} \\
(Y \mid Z = 1) & \sim \text{Poisson}(200)                              & \qquad \text{(target component)} \\
(Y \mid Z = 2) & \sim \text{Poisson}(70+\tau)                          & \qquad \text{(spillover component with shift)} \\
Y              & = (1-I) \cdot (Y \mid Z = 1) + I \cdot (Y \mid Z = 2) & \qquad \text{(mixture)}.
\end{aligned}
$$
The generative model for beads is an independent copy of the unshifted $Y \mid Z = 2$ at $\tau = 0$.

## Model Misspecification {-}

Generative model for real cells $Y$ of this experiment:
$$
\begin{aligned}
I              & \sim \text{Bernoulli}(0.1)                            & \qquad \text{(spillover indicator)} \\
Z              & = I + 1                                               & \qquad \text{(channel number)} \\
T              & \sim \text{Poisson}(200)                              & \qquad \text{(target)} \\
S              & \sim \text{Poisson}(70)                               & \qquad \text{(spillover)} \\
M              & \sim \text{Bernoulli}(\tau)                           & \qquad \text{(misspecification indicator)} \\
(Y \mid Z = 1) & = (1-M) \cdot T + M \cdot S                           & \qquad \text{(target mixture component)} \\
(Y \mid Z = 2) & = (1-M) \cdot S + M \cdot T                           & \qquad \text{(spillover mixture component)} \\
Y              & = (1-I) \cdot (Y \mid Z = 1) + I \cdot (Y \mid Z = 2) & \qquad \text{(mixture)}
\end{aligned}
$$
The generative model for beads is an independent copy of $Y \mid Z = 2$.

## Bimodal Spillover {-}

Generative model for real cells $Y$ of this experiment:
$$
\begin{aligned}
I               & \sim \text{Bernoulli}(0.1)                            & \qquad \text{(spillover indictor)} \\
Z               & = I + 1                                               & \qquad \text{(channel number)} \\
(Y \mid Z = 1)  & \sim \text{Poisson}(200)                              & \qquad \text{(target component)} \\
H               & \sim \text{Bernoulli}(\tau)                           & \qquad \text{(high count indicator)} \\
(S \mid H = 0)  & \sim \text{Poisson}(70)                               & \qquad \text{(low count component)} \\
(S \mid H = 1)  & \sim \text{Poisson}(330)                              & \qquad \text{(high count component)} \\
( Y \mid Z = 2) & = (1-H) \cdot (S \mid H = 0) + H \cdot (S \mid H = 1) & \qquad \text{(spillover component)} \\
Y               & = (1-I) \cdot (Y \mid Z = 1) + I \cdot (Y \mid Z = 2) & \qquad \text{(mixture)}
\end{aligned}
$$
The generative model for beads is an independent copy of $Y \mid Z = 2$.
