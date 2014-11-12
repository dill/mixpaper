---
title: PLOS revision checklist
author: David L Miller
---


# Comments to journal from authors

We than the AE and reviewers for their comments, which we have addressed line-by-line below. We do however believe that many of reviewer 1's comments did not indicate that the reviewer had sufficient understanding of the topic or had read the paper in detail. We would recommend that the journal not use this reviewer in the future.




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

* In the introduction, it would be clearer if the authors explain why they use finite mixtures functions rather than continuous mixture functions. It would be worthwhile to explain the differences between finite and continuous mixture models.
  - *We go on to explain the continuous mixture idea in the discussion. Including this information in the introduction and then going on to not use continuous mixtures is potentially misleading for readers. We have added an explanation of why we do not use continuous mixtures in the discussion and clearly state in the introduction that we are using finite mixtures.*
* Throughout the manuscript, it was difficult to me to discern the motivation of the simulation study. Some additional explanations that describe explicitly the aim of the simulation study are needed in the introduction and the section that describes the simulation study.
  - *From the "Examples: Simulated data section": "Extensive simulations were carried out to investigate performance (in terms of the accuracy of estimation of $P_a$) when the true detection function model is not known to the estimation procedure." Within that section each simulation scenario is also justified.*
* Although, the main motivation of the paper is based on the need of monotone and non-increasing functions to model distributions of probabilities of detection. I would have like to read more arguments about the importance of introducing mixture models in distance sampling framework. Indeed, I am unclear about how far this new class of functions can be applied. Is it specific to given organisms, to a particular distribution of probability of detection or is it useful for any organisms for which distance sampling data have been collected? It would be worthwhile that the authors clarify these points.
  - *Monotone detection functions are the main reason that these models have been introduced. As with all detection functions, these can be applied to any species.*
* In the introduction, the authors evoke the possibility of using mixture models which account for highly heterogeneous detection probability. Why do they not propose such models in this paper? It would be interesting that the authors propose a simulation case where the finite mixture models are fitted to probabilities of detection that are highly heterogeneous. It could be helpful to assess how the finite mixture models estimate the probability of detection in such extreme case.
  - *Indeed, we do propose simulations in which one component is a spike. A4/B4 and C2 provide good examples of this.*
* If the authors do not plan to propose combination of continuous and finite mixture models in the future version of this manuscript, they should move the part of the introduction (page 3 from "In addition, mixture models" ... to "to an appealing conceptual explanation for underlying data") in the discussion.
  - *This comment does not make sense. Combination of continuous and finite mixtures is not the only way to deal with heterogeneity. As we show in simulations, finite mixtures can deal with this situation.*


### Specific comments

#### Methods

* "K-vector of the associated covariates", more explanations about what is a "K-vector" would be helpful.
  - *this is a standard way to denote a vector of length $K$*
* Why the authors assume that each mixture component has a different scale, Is not just an intrinsic property of finite mixture models?
  - *Think this is down to over-interpretation on the part of the reviewer, clarified.*
* I wonder if there is not a problem of notation in equations 2 and 3. The scale parameter $\sigma_j$ in equation 2 becomes $\sigma_{ij}$ in equation 3. In equation 3, a subscript $i$ is added, as a consequence the link between equation 2 and 3 becomes unclear. If it is not a problem of notation, more explanations about "K-vector" are needed.
  - *The equations are "linked" as in Marques and Buckland (2003). Equation 2 is not index by $i$, so the additional subscipt would not make sense. As indicated in the text, $i$ indexes the observations, so $\sigma_{ij}$ is a shorthand for the full exponential form for the scale parameter for a given observation and it's associated covariates.*
* What is the truncation distance? More explanations on the truncation distance are needed, especially on how the truncation distance is chosen.
  - *Truncation distance is covered in the first paragraph of the introduction, added minor comment at this point. Choice of truncation distance is out of scope for this article, but is covered thoroughly throughout the distance sampling literature (e.g. Buckland et al (2001), now referenced at that point).*
* I am also wondering how measurement errors of covariates might affect the estimation of the probability of detection. Is there any specific kind of covariates that can be introduced in the finite mixture models? If a covariate is measured on a transect while it affects the individuals at the place where they are observed, how bias it the estimation of the probability of detection function? Is it a large source of bias to use covariates that are not measured at the place where the individuals are observed. Some words or some references could be added on how the localisation of the measurement of a covariate can affect the estimation of the probability of detection.
  - *Measurement error for covariates is out of the scope of this article, though we note that covariates are effectively nuisance variables ther we use in order to explain heterogeneity in detectability, we don't perform inference using them. Measurement error for distances has been investigated in the literature (e.g. Borchers et al, 2010); we have added this into the discussion.*
  - *As in Marques and Buckland (2003), there are no restrictions on covariates that are to be included in the detection function.*
  - *Since we are interested in covariates that affect detection, we usually don't have to worry about the mismatch between the covariate at the observer's location and the observations location. Many commonly collected covariates are not afftected by proximity (e.g. sex, behaviour or group size), though their correct identification might be. Others such at sea state or visibility will affect detectability where the observer **is**, not where the object is observed.*

#### Estimating population size

* The authors give the equations to estimate population size while they do not provide any estimate of population size either in the simulated study or for the case studies. It would also be interesting to have the estimates of the population size (N) for the case studies.
  - *Abundance estimates were omitted from the case studies and simulation results as we did not want to further complicate the presentation of the results by including a further metric for readers to consider. For the simulations, the abundance estimates in the covered area are simply functions of the probabilities presented. For the case studies, data come from relatively complex survey designs, so meaningful abundance estimates would require significantly more exposition. For these reasons we choose to ommit estimated of abundance.*

#### Simulated data

* The authors propose to use half-normal distribution, exponential power series, two hazard-rate function, but it remains unclear throughout the manuscript why they use such functions (except that they are non-increasing and monotone). Does it possible to relate the different functions to "biological reality"? if yes or even no, it should be better explained in this part of the Ms.
  * *The half-normal and hazard-rate functions are used throughout the distance sampling literature as they fullfil the criteria listed in Buckland et al (2001) for detection functions: that functions are flexible but use a relatively small number of parameters, while conforming to the shape we expect from examining empirical data. (Section 2.3 of Buckland et al, 2001 covers this in depth.*
* In Group C (page 6): why the authors want to add "complication" and that "one of the components has a larger scale parameter relative to truncation distance"?
  - *The "complication" is added to make the simulation more complex and hence test the model in a more difficult situation. In this case, when the scale parameter large relative to the truncation, the exact size of the scale is hard to estimate as the effect of a scale parameter of 10 and 100 is not distinguishable.*

#### Case studies

* In the abstract, the authors state "We also re-analyze four previously problematic real-world case studies" While in the main text, in the section where they describe the case studies (page 6), they do not clearly explain why the case studies are problematic. It is only in the legend of the Figure 1 that it is explained clearly that the detection functions are non-monotone for two examples. If the authors add explanations on the problems encountered with each of the case studies (in the section where they describe the case studies), it will clarify the motivation of the choice of the case studies.
  - *Information on why the data sets were included is gvein in the "Case studies" section on page 6/7*.

#### Discussion

* P 9. "we note that is better to avoid collecting such data in the first place" Why it is better to avoid collecting such data? How can it be avoided to collect such data? The reference that is given, there, is a book. It would be worth at least to give the pages of the book that deal with this problem.
  - *Explanation of why spiked data is problematic is given on page 3. Added pages to reference. Spiked data indicate that observers have been spending too much time looking on or very near the trackline to the detriment of objects further away.*
* P 10. From "A half-normal detection function is relatively inflexible so uncertainty it low" to "Mixtures of half-normals lie somewhere in between these two options" Why the inflexibility of the half-normal function yields low uncertainty? Some additional explanations are needed there.
  - *Corrected "it"/"is" typo. Additional explanation in terms of the scale parameter have been added.*
* In the fourth paragraph of page 10, (i.e. from "In simulation we observed"), the sentence "2-point mixtures were generally chosen by AIC as good models (However these models were useful". I am unclear about which models are useful, also why they are useful.
  - *Additional text has been added here to clarify.*
* At the end of the discussion, the authors suggest other component functions that could be useful. Why they do not try to use it in this Ms? and in which case it would be more useful than half-normal function?
  * *The other obvious candidate for mixture component would be the hazard-rate function. This was not pursued in this manuscript as we had believed we had already introduced a rather large amount of new material and presented a lot of simulation and case study results. We are also unsure of the additional utility of the hazard-rate function and exactly what the formulation would consist of. For example would a common shape parameter be appropriate between components, or would separate scales be more useful. We also note that a hazard-rate mixture in the latter formulation would be adding 3 parameters per mixture component, so a 2-point mixture of hazard-rate functions would use 5 parameters without any covariates in the model. In K+A approaches we rarely select models with >5 parameters, so hazard-rate models seem uncompetative from an AIC model selection perspective.*
* In which case a combination of both finite and continuous mixtures functions could be used? Why the authors do not use it in this paper? Some explanations about the difference between continuous and finite mixtures functions are needed in the Ms. Why does it seem important to the authors to echoes the work on mixture models that is done in capture-recapture framework?
  - *Continuous mixtures present a much more complex and time consuming optimisation problem, a mixture of finite and continuous mixtures is even more complex from an optimisation perspective. In keeping with the comments above, we also believed that there was already a rather large amount of new material in the paper already without adding more complicated modelling options (and their associated simulations etc).*

#### Table and Figures

 * Table 1: How the different models are stored in the table is unclear. It would be clearer to me, if, for each case study, the models were ranked by increasing values of $\Delta$AIC.
  - *The first item in each case study group is the "original" model fitted using K+A methods in the previous analysis of the data. We have now sorted the other models in each case study group according to their $\Delta$AIC.*
 * Figure 1: Why the shape of the functions of A2, A4, C1 and D1 are so similar? Why there are 3 dashed lines in C1 and C2? Why there are no dashed lines in E1?
  - *Assuming the reviewer is referring to figure 2...*
  - *We conend that A4 is not similar to the other three functions, as it has a rather more pronounced spike than the others. Functions A2, C1 and D1 represent rather different situations -- A2 a simple 2-point mixture, C1 a 3-point mixture that looks like it could be well approximated by a 2-point mixture and D1 is a covariate model, so the detection function give is presented at the median levels of the covariates and actually represents the median effect.*
  - *As in the previous figures, the dashed lines in C1 and C2 indicate the components. We have clarified in the caption.*
  - *E1 has no dashed lines as it does not have any componets, it is just the exponential power series model.*
 * Figure 3: In the second panels of scenario D and scenario E, why the values of the proportion of AIC best models that were of the same form as the model that the data was simulated from, are not given?
  - *Models of the same form as the simulated model were not fitted to the data in scentarios E1 and E2. As stated in the main text, the second set of panels for scenario D do not include the covariates in order to investigate whether mixture components could take the place of covariates -- so no models that were of the same form as the data generation process were fitted.*

## Reviewer 2

  * page 8, line 10 from bottom: slightly
    - *Fixed.*
  * page 9, line 4-5: when P_a is lower, wouldn't that cause an underestimate of N, rather than an overestimate as stated ? Whatever the case, you may rewrite this bit to avoid confusion.
    - *We agree that this statement was poorly worded and incorrect, we have changed the text accordingly.*
  * page 9, line 7 from bottom: perhaps not sooo surprising if one considers that the K+A approach dates back more than 20 years. Other areas of statistics have come a long way since. One might say it is surprising that the old approach has 'survived' for so long.
    - *We have removed the middle part of this "that the K+A formulation is designed to produce parsimonious and realistic fits," as we do still think it's surprising that models were selected even though they had more parameters than K+A models.*
  * page 10, line 7: say variation in what
    - *Clarified that we mean CV of the average probability of detection.*
  * same page, middle of page, "however these models": make clear which one you mean
    - *Can't find this reference.*
  * just below that: I think that in the finite-mixture work of Pledger (e.g., 2000) she also finds that mixtures with more than 2 or 3 support points are rarely supported by AIC. Or, I believe they often won't converge when one tries to fit them.
    - *Indeed, Pledger (2000) finds that usually only 2-point mixtures are selected. We have added comments to this effect to the text.*
  * One main (but really small) criticism I have is that I found neither the tables nor the figures very nicely produced. For instance, in Table 1, I would suggest the following changes. give more informative table legend, where the meanuing of each column heading is explained (even if is in the text already, I find it good when a table can be understood without referring to the main text); make a different separation for different models fit to the same data set and those fitted to different data sets; perhaps print the AIC-minimal model in bold face; print 2 digits overall. Fig. 3 was way too crowded for my taste, perhaps I distribute the results onto more than one table; then, I would show the axes of each diagram
    - *Each column heading is now explained in the caption.*
    - *We wonder what is possible from the journal's point of view -- would it be possible to have species highlighted in alternating colours?*
    - *As per reviewer 1's comments, we have sorted the models according to AIC. We have also expanded all model's AICs to 2 decimal places.*
    - *Though we agree that Figure 3 contains a lot of information, we do think that it is useful for readers to be able to see the proportion of top models that were of the same form as the model used to generate the data.*
  * Appendix S1, line 4 in main text: Shannon; the axis labels in the Figure 1 in same App are too small
    * *D. F. Shanno's surname is "Shanno", see e.g. Nocedal & Wright (2000) section 6.1.*
    - *We have increased the size of the figure and text in Appendix S1, figure 1.*
  * At the start of Appendix S2, I would remind the reader what this is all about, e.g., "We conducted a simulation with ....". Again, I think you can make the appendix understandable by itself at little cost, e.g., by adding 1-2 sentences not more. Same place, Table 1: explain what all the column headings mean.
    - *We have added explanation to the start of Appendix S2.*
    - *Added explanation of column headings.*


# References

  * Borchers, D., Marques, T., Gunnlaugsson, T., & Jupp, P. (2010). Estimating Distance Sampling Detection Functions When Distances Are Measured With Errors, 15(3), 346–361. doi:10.1007/s13253-010-0021-y
  * Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L., Borchers, D. L., & Thomas, L. (2001). Introduction to Distance Sampling. Oxford University Press.
  * Pledger, S. “Unified Maximum Likelihood Estimates for Closed Capture-Recapture Models Using Mixtures.” Biometrics 56, no. 2 (2000): 434–42.
  * Nocedal, J., Wright, S. J. (2000). Numerical Optimization. Springer.
