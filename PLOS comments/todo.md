---
title: PLOS revision checklist
author: David L Miller
---

# Journal requirements

1. [ ] If your submission was prepared in LaTex, please submit your manuscript file in PDF format and attach your .tex file as “other.”
2. [ ] We note that you stated “data are available upon request” at submission. Could you please confirm that all data underlying the findings in your study are freely available in the manuscript, supplemental files, or in a public repository? If this is not the case, and your data are available upon request because of an ethical or legal restriction or because you obtained data from a third party, please include the following in your revised cover letter:
  a. [ ] The reason why your data cannot be made available in the manuscript, the supplemental files, or a public repository;
  b. [ ] The name(s) of the individual(s) that readers may contact to request the data;
  We will make changes to your data availability statement on your behalf, based on the information you provide. For more information about our data policy and acceptable reasons for not making your data fully available, please refer to: http://www.plosone.org/static/policies#sharing
3. [ ] We also note that you have stated that you will provide repository information for your data at acceptance. Should your manuscript be accepted for publication, we will hold your manuscript until you get in touch with us with the accession numbers or DOIs necessary to access your data. If you wish to make changes to your data availability statement, please describe these changes in your cover letter and we will make them on your behalf.

# From AE

I agree with most of their (minor) remarks, which, if you can take them into account, could help clarify the manuscript. **Please notably consider the comments of referee 2 regarding the figures and tables.**

# Reviewer comments

## Reviewer 1

### General comments

 * [ ] In the introduction, it would be clearer if the authors explain why they use finite mixtures functions rather than continuous mixture functions. It would be worthwhile to explain the differences between finite and continuous mixture models.
 * [ ] Throughout the manuscript, it was difficult to me to discern the motivation of the simulation study. Some additional explanations that describe explicitly the aim of the simulation study are needed in the introduction and the section that describes the simulation study.
 * [ ] Although, the main motivation of the paper is based on the need of monotone and non-increasing functions to model distributions of probabilities of detection. I would have like to read more arguments about the importance of introducing mixture models in distance sampling framework. Indeed, I am unclear about how far this new class of functions can be applied. Is it specific to given organisms, to a particular distribution of probability of detection or is it useful for any organisms for which distance sampling data have been collected? It would be worthwhile that the authors clarify these points.
 * [ ] In the introduction, the authors evoke the possibility of using mixture models which account for highly heterogeneous detection probability. Why do they not propose such models in this paper? It would be interesting that the authors propose a simulation case where the finite mixture models are fitted to probabilities of detection that are highly heterogeneous. It could be helpful to assess how the finite mixture models estimate the probability of detection in such extreme case.
 * [ ] If the authors do not plan to propose combination of continuous and finite mixture models in the future version of this manuscript, they should move the part of the introduction (page 3 from "In addition, mixture models" ... to "to an appealing conceptual explanation for underlying data") in the discussion.

### Specific comments

#### Methods

 * [ ] "K-vector of the associated covariates", more explanations about what is a "K-vector" would be helpful.

 * [ ] Why the authors assume that each mixture component has a different scale, Is not just an intrinsic property of finite mixture models?

 * [ ] I wonder if there is not a problem of notation in equations 2 and 3. The scale parameter $\sigma_j$ in equation 2 becomes $\sigma_{ij}$ in equation 3. In equation 3, a subscript $i$ is added, as a consequence the link between equation 2 and 3 becomes unclear. If it is not a problem of notation, more explanations about "K-vector" are needed.

 * [ ] What is the truncation distance? More explanations on the truncation distance are needed, especially on how the truncation distance is chosen.

 * [ ] I am also wondering how measurement errors of covariates might affect the estimation of the probability of detection. Is there any specific kind of covariates that can be introduced in the finite mixture models? If a covariate is measured on a transect while it affects the individuals at the place where they are observed, how bias it the estimation of the probability of detection function? Is it a large source of bias to use covariates that are not measured at the place where the individuals are observed. Some words or some references could be added on how the localisation of the measurement of a covariate can affect the estimation of the probability of detection.

#### Estimating population size

 * [ ] The authors give the equations to estimate population size while they do not provide any estimate of population size either in the simulated study or for the case studies. It would also be interesting to have the estimates of the population size (N) for the case studies.


#### Simulated data

 * [ ] The authors propose to use half-normal distribution, exponential power series, two hazard-rate function, but it remains unclear throughout the manuscript why they use such functions (except that they are non-increasing and monotone). Does it possible to relate the different functions to "biological reality"? if yes or even no, it should be better explained in this part of the Ms.

 * [ ] In Group C (page 6): why the authors want to add "complication" and that "one of the components has a larger scale parameter relative to truncation distance"?

#### Case studies

 * [ ] In the abstract, the authors state "We also re-analyze four previously problematic real-world case studies" While in the main text, in the section where they describe the case studies (page 6), they do not clearly explain why the case studies are problematic. It is only in the legend of the Figure 1 that it is explained clearly that the detection functions are non-monotone for two examples. If the authors add explanations on the problems encountered with each of the case studies (in the section where they describe the case studies), it will clarify the motivation of the choice of the case studies.


#### Discussion

 * [ ] P 9. "we note that is better to avoid collecting such data in the first place" Why it is better to avoid collecting such data? How can it be avoided to collect such data? The reference that is given, there, is a book. It would be worth at least to give the pages of the book that deal with this problem.

 * [ ] P 10. From "A half-normal detection function is relatively inflexible so uncertainty it low" to "Mixtures of half-normals lie somewhere in between these two options" Why the inflexibility of the half-normal function yields low uncertainty? Some additional explanations are needed there.

 * [ ] In the fourth paragraph of page 10, (i.e. from "In simulation we observed"), the sentence "2-point mixtures were generally chosen by AIC as good models (However these models were useful". I am unclear about which models are useful, also why they are useful.

 * [ ] At the end of the discussion, the authors suggest other component functions that could be useful. Why they do not try to use it in this Ms? and in which case it would be more useful than half-normal function?

 * [ ] In which case a combination of both finite and continuous mixtures functions could be used? Why the authors do not use it in this paper? Some explanations about the difference between continuous and finite mixtures functions are needed in the Ms. Why does it seem important to the authors to echoes the work on mixture models that is done in capture-recapture framework?

#### Table and Figures

 * [ ] Table 1: How the different models are stored in the table is unclear. It would be clearer to me, if, for each case study, the models were ranked by increasing values of $\Delta$AIC.

 * [ ] Figure 1: Why the shape of the functions of A2, A4, C1 and D1 are so similar? Why there are 3 dashed lines in C1 and C2? Why there are no dashed lines in E1?

 * [ ] Figure 3: In the second panels of scenario D and scenario E, why the values of the proportion of AIC best models that were of the same form as the model that the data was simulated from, are not given?

## Reviewer 2

 * [ ] page 8, line 10 from bottom: slightly
 * [ ] page 9, line 4-5: when P_a is lower, wouldn't that cause an underestimate of N, rather than an overestimate as stated ? Whatever the case, you may rewrite this bit to avoid confusion.
 * [ ] page 9, line 7 from bottom: perhaps not sooo surprising if one considers that the K+A approach dates back more than 20 years. Other areas of statistics have come a long way since. One might say it is surprising that the old approach has 'survived' for so long.
 * [ ] page 10, line 7: say variation in what
 * [ ] same page, middle of page, "however these models": make clear which one you mean
 * [ ] just below that: I think that in the finite-mixture work of Pledger (e.g., 2000) she also finds that mixtures with more than 2 or 3 support points are rarely supported by AIC. Or, I believe they often won't converge when one tries to fit them.
 * [ ] One main (but really small) criticism I have is that I found neither the tables nor the figures very nicely produced. For instance, in Table 1, I would suggest the following changes. give more informative table legend, where the meanuing of each column heading is explained (even if is in the text already, I find it good when a table can be understood without referring to the main text); make a different separation for different models fit to the same data set and those fitted to different data sets; perhaps print the AIC-minimal model in bold face; print 2 digits overall. Fig. 3 was way too crowded for my taste, perhaps I distribute the results onto more than one table; then, I would show the axes of each diagram
 * [ ] Appendix S1, line 4 in main text: Shannon; the axis labels in the Figure 1 in same App are too small
 * [ ] At the start of Appendix S2, I would remind the reader what this is all about, e.g., "We conducted a simulation with ....". Again, I think you can make the appendix understandable by itself at little cost, e.g., by adding 1-2 sentences not more. Same place, Table 1: explain what all the column headings mean.


