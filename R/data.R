#' Custom color palettes
#'
#' @description
#' A collection of hand-picked, coordinated color palettes.
#' These palettes are specifically designed for high-contrast, publication-ready
#' data visualizations. Recently expanded to 7 colors per palette to accommodate
#' datasets with more categorical levels.
#'
#' @format A named list containing 8 vectors. Each vector contains 7 named hex color codes:
#' \describe{
#'   \item{maritime}{Dark blues contrasted with warm terracotta, rounded out by kelp green and soft sand.}
#'   \item{botanical}{Earthy greens and deep berry tones, balanced with soft rose and sage.}
#'   \item{nordic}{Teal (petrol) and crimson, grounded by deep pine and light linen.}
#'   \item{azure_coral}{Bright azure blue and fresh coral red, with deep navy and soft peach.}
#'   \item{eucalyptus}{Eucalyptus green and raspberry red, accented by dark fir and light cream.}
#'   \item{aegean}{Modern Mediterranean look with ocean blue and sand tones, deepened by dark sea and cloud white.}
#'   \item{oasis}{Strong teal and grapefruit contrast, extended with deep midnight blue and soft peach.}
#'   \item{iris_copper}{Elegant cornflower blue and copper, enriched with aubergine and warm sand.}
#' }
"custom_palettes"

#' Simulated Thesis Writing Progress Data
#'
#' @description
#' A simulated dataset in long format reflecting the writing productivity
#' of individuals. It documents the total number of words written at two
#' consecutive time points: before and after a designated writing period
#' or intervention.
#'
#' @details
#' **Academic Context & Meaningful Progress:**
#' In writing research, tracking raw word counts over time provides insight
#' into writer productivity and fluency. While any positive change indicates
#' progress, evaluating the efficacy of a writing intervention often requires
#' defining a specific threshold. For instance, an individual improvement of
#' at least 50 words might be considered a meaningful short-term gain in a
#' micro-writing session context, distinguishing highly productive phases
#' from natural, minor fluctuations.
#'
#' **Background of the Simulation:**
#' The data simulates a large-scale pre-post study with 10,000 participants,
#' resulting in 20,000 total observations. The simulation parameters mirror
#' highly variable individual writing habits.
#'
#' The data generating process was defined as follows:
#' \itemize{
#'   \item Baseline ("Before") values are normally distributed around a mean of 600 words (SD = 150).
#'   \item The true individual change (treatment effect) was simulated with a mean increase of 20 words (SD = 50), reflecting that while average growth is modest, individual variations are substantial.
#'   \item "After" values represent the sum of the baseline values and the individual effects.
#' }
#'
#' @format A tibble with 20,000 rows and 3 variables:
#' \describe{
#'   \item{id}{Unique participant identification number (1 to 10000).}
#'   \item{time}{Measurement time point: \code{"Before"} or \code{"After"}.}
#'   \item{words_written}{Total number of words written recorded at the specific time point.}
#' }
"thesis"

#' Simulated Frozen Shoulder Range of Motion Data
#'
#' @description
#' A simulated dataset in long format reflecting the treatment outcomes of
#' patients with "frozen shoulder" (adhesive capsulitis). It documents the
#' maximum active range of motion (ROM) of the shoulder in degrees at two
#' time points: before (pre) and after (post) a physical therapy intervention.
#'
#' @details
#' **Clinical Context & Minimal Important Difference (MID):**
#' In clinical practice, a purely statistical improvement in range of motion
#' is often insufficient to alleviate a patient's burden. To experience a
#' noticeable relief in everyday activities (e.g., brushing hair, dressing),
#' a specific threshold must be crossed. For this dataset, a **Minimal
#' Important Difference (MID) of 15 degrees** is defined. An individual
#' improvement of at least 15 degrees is therefore considered a
#' clinically meaningful treatment success.
#'
#' **Reliability & Measurement Error:**
#' The reliability of the measurement instrument (e.g., a universal goniometer)
#' is a crucial parameter for Bayesian individual change models. For this
#' dataset, a **test-retest reliability (ICC) of 0.90** is assumed, reflecting
#' a standardized intra-rater scenario.
#'
#' This reliability is used to calculate the Standard Error of the Difference
#' (\eqn{S_{diff}}), which serves as the measurement error (sigma) in the
#' Bayesian engine:
#' \deqn{S_{diff} = sd(pre\_values) \times \sqrt{2 \times (1 - r_{xx})}}
#'
#' **Background of the Simulation:**
#' The data simulates a typical pre-post study (n = 50, observations = 100).
#' Baseline values are normally distributed around a mean of 90 degrees (SD = 15).
#' True change per patient was simulated with a mean of 20 degrees (SD = 12).
#'
#' @format A tibble with 100 rows and 3 variables:
#' \describe{
#'   \item{id}{Unique patient identification number (1 to 50).}
#'   \item{time}{Measurement time point: "Pre" (baseline) or "Post" (follow-up).}
#'   \item{rotation}{Range of motion (ROM) in degrees.}
#' }
#'
#' @source Simulated for demonstration purposes of Clinical Significance Analyses.
"shoulder"


#' Simulated Emotional Exhaustion Inventory (EEI) Data
#'
#' @description
#' A simulated longitudinal dataset representing Emotional Exhaustion Inventory (EEI)
#' scores in a clinical trial comparing a Waitlist group to an Intervention group.
#' The data includes baseline measurements and follow-ups over varying time points.
#'
#' @details
#'
#' For clinical decision making, a \strong{Minimal Important Difference (MID)}
#' of 10 points is assumed. This leads to a \strong{Region of Practical Equivalence (ROPE)}
#' of \code{[-10, 10]}.
#' \itemize{
#'   \item \strong{Clinical Improvement:} A reduction in EEI scores of > 10 points.
#'   \item \strong{Clinical Deterioration:} An increase in EEI scores of > 10 points.
#' }
#'
#' @format A data frame (tibble) with 261 rows and 4 variables:
#' \describe{
#'   \item{id}{Unique identifier for each participant (1 to 60).}
#'   \item{group}{Factor with two levels: \code{"Waitlist"} and \code{"Intervention"}.}
#'   \item{time_weeks}{Time since baseline in weeks. Intervention group follow-ups
#'     vary (2-6 weeks), while Waitlist group follows a fixed schedule (0, 4, 8, 12 weeks).}
#'   \item{eei}{Emotional Exhaustion Inventory score (0 to 100). Higher scores
#'     indicate higher levels of burnout. Lower scores are benefitial.}
#' }
#'
#' @usage data(burnout)
#'
#' @source Simulated data generated for clinical certainty analysis.
#'
#' @examples
#' data(burnout)
"burnout"


#' Simulated Morning Energy and Espresso Data
#'
#' @description
#' A clean, simulated longitudinal dataset designed to demonstrate hierarchical Bayesian models,
#' "smart shrinkage", and clinical certainty. It tracks the morning energy levels
#' of 50 patients over 4 weeks, with an intervention starting in week 2.
#'
#' @details
#' The dataset was generated to illustrate how contextual knowledge (e.g., espresso
#' consumption) can help accurately estimate the true effect of an intervention when
#' dealing with highly noisy measurements.
#'
#' While the underlying data generation process includes individual baseline deviations
#' and varying intervention effects, this dataset only exposes the observed variables.
#' A naive model might be easily confused by the measurement noise (SD = 4), whereas a
#' hierarchical model can reveal the underlying structure through shrinkage.
#'
#' For clinical decision making, a \strong{Minimal Important Difference (MID)}
#' of 5 points is assumed, resulting in a \strong{Region of Practical Equivalence (ROPE)}
#' of \code{[-5, 5]}.
#' \itemize{
#'   \item \strong{Clinical Success:} An increase in energy of > 5 points.
#'     Patients consuming 0 to 2 espressos exceed this threshold in the true data generation process,
#'     while patients consuming 3 or 4 espressos experience no clinically relevant improvement in reality.
#' }
#'
#' @format A data frame (tibble) with 200 rows (50 patients x 4 weeks) and 6 variables:
#' \describe{
#'   \item{id}{Unique identifier for the patient (1 to 50).}
#'   \item{age}{Age of the patient in years (approx. 20 to 60).}
#'   \item{espresso}{Habitual daily espresso consumption (0 to 4).}
#'   \item{week}{Week of measurement (0, 1, 2, 3).}
#'   \item{is_post}{Binary indicator for the intervention phase (0 = weeks 0-1, 1 = weeks 2-3).}
#'   \item{morning_energy}{Observed (measured) energy score including realistic measurement noise.}
#' }
#'
#' @usage data(espresso)
#'
#' @source Simulated dataset for methodological demonstrations.
#'
#' @examples
#' data(espresso)
"espresso"


#' Simulated Single-Case Data for Misophonia Treatment
#'
#' @description
#' A simulated Single-Case Experimental Design (SCED) dataset tracking the daily
#' misophonia trigger severity of a single patient over 30 days. The data is divided
#' into a baseline phase and an intervention phase.
#'
#' @details
#' This dataset illustrates a classic A-B phase design for a rare or highly
#' individualized psychological condition. The patient suffers from misophonia
#' (extreme emotional reactivity to specific trigger sounds).
#'
#' The underlying data generation process assumes:
#' \itemize{
#'   \item \strong{Baseline (Days 1-10):} Without treatment, the patient experiences
#'     daily sensitization, leading to a steady increase in severity scores (slope > 0).
#'   \item \strong{Therapy (Days 11-30):} The introduction of a specific habituation
#'     intervention reverses the trend, causing a sharp daily decrease in severity (slope < 0).
#'   \item \strong{Autoregressive Noise:} Daily measurements are not independent. An AR(1)
#'     process simulates how a particularly bad day negatively impacts the mood of the following day.
#' }
#'
#' Severity scores are bounded between 0 (no distress) and 60 (extreme distress).
#'
#' @format A data frame (tibble) with 30 rows and 3 variables:
#' \describe{
#'   \item{day}{The day of the study (1 to 30).}
#'   \item{phase}{The study phase. Either \code{"Baseline"} (days 1-10) or \code{"Therapy"} (days 11-30).}
#'   \item{severity_score}{Daily distress level upon encountering the primary audio trigger (0 to 60).}
#' }
#'
#' @usage data(misophonia)
#'
#' @source Simulated dataset for methodological demonstrations of SCEDs.
#'
#' @examples
#' data(misophonia)
"misophonia"


#' Simulated Single-Case A-B-A-B Design for Erythromelalgia
#'
#' @description
#' A simulated Single-Case Experimental Design (SCED) dataset tracking the daily
#' pain levels of a single patient over 40 days using an A-B-A-B withdrawal design.
#'
#' @details
#' This dataset illustrates a classic A-B-A-B reversal design applied to a rare
#' and severe condition. The patient suffers from treatment-resistant Erythromelalgia
#' (often called "Man-on-Fire Syndrome"), caused by a mutation in the Nav1.7 sodium channel.
#' Due to the rarity of the condition and the failure of standard treatments, an
#' individual compassionate use trial (N=1) is conducted using an experimental
#' Nav1.7 channel blocker (a topical cone-snail peptide analog).
#'
#' The underlying data generation process highlights specific pharmacological dynamics:
#' \itemize{
#'   \item \strong{Phase A1 (Baseline, Days 1-10):} Severe, unmanaged burning pain (Mean = 80).
#'   \item \strong{Phase B1 (Intervention, Days 11-20):} The topical peptide blocks the channels,
#'     leading to immediate and massive pain relief (Mean = 35).
#'   \item \strong{Phase A2 (Withdrawal, Days 21-30):} The gel is removed to prove causality.
#'     Due to receptor upregulation during the blockade, a "rebound effect" occurs,
#'     making the pain slightly worse than the initial baseline (Mean = 85).
#'   \item \strong{Phase B2 (Re-Intervention, Days 31-40):} The gel is reapplied,
#'     promptly restoring pain relief and proving the causal efficacy of the drug (Mean = 30).
#'   \item \strong{Autoregressive Noise:} An AR(1) process simulates carry-over effects,
#'     meaning a particularly bad pain day slightly elevates the baseline for the next day.
#' }
#'
#' @format A data frame (tibble) with 40 rows and 3 variables:
#' \describe{
#'   \item{day}{The day of the study (1 to 40).}
#'   \item{phase}{The study phase: \code{"A1"} (Baseline), \code{"B1"} (Intervention), \code{"A2"} (Withdrawal), or \code{"B2"} (Re-Intervention).}
#'   \item{pain_score}{Daily reported pain severity bounded between 0 (no pain) and 100 (worst imaginable pain).}
#' }
#'
#' @usage data(man_on_fire)
#'
#' @source Simulated dataset for methodological demonstrations of A-B-A-B designs.
#'
#' @examples
#' data(man_on_fire)
"man_on_fire"


#' Simulated Return-to-Work Data (Sickness Absence)
#'
#' @description
#' A simulated longitudinal dataset tracking the number of monthly sick days
#' for 50 employees over a period of 6 months. This dataset is designed to
#' demonstrate Negative Binomial regression and hierarchical modeling of
#' count data with overdispersion.
#'
#' @details
#' The data generation process follows a log-linear trend where the expected
#' number of sick days decreases over time.
#'
#' Key characteristics of the simulation:
#' \itemize{
#'   \item \strong{Baseline:} Employees start with an average of approximately 12
#'     sick days in Month 1 (\code{exp(log(12))}).
#'   \item \strong{Recovery Rate:} On average, the number of sick days drops by
#'     20% each month (\code{exp(log(0.80))}), representing a successful
#'     return-to-work intervention or natural recovery.
#'   \item \strong{Individual Variability:} Random intercepts and random slopes
#'     account for the fact that some employees start with higher absence levels
#'     and some recover faster than others.
#'   \item \strong{Overdispersion:} Actual sick days are generated using a
#'     Negative Binomial distribution (\code{size = 100}) to allow for more
#'     variance than a standard Poisson model would permit.
#' }
#'
#' The \code{sick_days} variable is capped at 32 days per month to maintain
#' biological plausibility.
#'
#' @format A data frame (tibble) with 300 rows (50 employees x 6 months) and 4 variables:
#' \describe{
#'   \item{id}{Unique identifier for the employee (1 to 50).}
#'   \item{month}{Month of observation (1 to 6).}
#'   \item{time}{Zero-centered time variable (\code{month - 1}), where 0 represents the first month.}
#'   \item{sick_days}{The number of sick days recorded in that month (count data, 0-32).}
#' }
#'
#' @usage data(return_to_work)
#'
#' @source Simulated dataset for demonstrating Count Data GLMMs.
#'
#' @examples
#' data(return_to_work)
"return_to_work"


#' Simulated Multivariate Athlete Performance Data
#'
#' @description
#' A simulated longitudinal dataset of 60 athletes measuring physical strength
#' and performance anxiety before and after a specialized training camp.
#' It is specifically designed to demonstrate multivariate hierarchical models
#' and the recovery of correlated random effects.
#'
#' @details
#' This dataset illustrates a common challenge in sports psychology and biostatistics:
#' modeling two distinct but interacting outcomes simultaneously. The training camp
#' aims to increase physical strength while decreasing performance anxiety.
#'
#' Key features of the simulation:
#' \itemize{
#'   \item \strong{Correlated True Effects:} At the latent level, there is a strong
#'     negative correlation (\eqn{r = -0.75}) between an athlete's treatment response
#'     for anxiety and strength. Athletes who experience a larger drop in anxiety
#'     tend to see a greater gain in strength.
#'   \item \strong{Repeated Measurements:} To accurately recover this true correlation,
#'     daily measurement noise must be separated from true individual differences.
#'     Therefore, each athlete is measured 3 times before the camp and 3 times after.
#'   \item \strong{Measurement Noise:} The observed scores fluctuate around the
#'     true latent scores due to realistic daily variance (SD = 8 for strength,
#'     SD = 4 for anxiety).
#' }
#'
#' A naive correlation of simple difference scores will underestimate the true
#' relationship due to measurement error (attenuation). A multivariate mixed
#' model uses the replicate measurements to partial out the noise and recover
#' the true correlation of the random slopes.
#'
#' @format A data frame (tibble) with 360 rows (60 athletes x 2 phases x 3 replicates) and 5 variables:
#' \describe{
#'   \item{id}{Unique identifier for the athlete (1 to 60).}
#'   \item{phase}{Binary indicator for the phase (0 = Pre-camp, 1 = Post-camp).}
#'   \item{measurement}{Index of the repeated measurement within the phase (1, 2, or 3).}
#'   \item{strength}{Observed physical strength score, including measurement noise.}
#'   \item{anxiety}{Observed performance anxiety score, including measurement noise.}
#' }
#'
#' @usage data(athletes)
#'
#' @source Simulated dataset for methodological demonstrations of multivariate models.
#'
#' @examples
#' data(athletes)
"athletes"


#' Simulated Post-Cardiac Rehabilitation Data (Missing At Random)
#'
#' @description
#' A simulated longitudinal dataset tracking the physical activity of 90 patients
#' following a myocardial infarction (heart attack) over a 12-wave rehabilitation
#' program. It is specifically designed to demonstrate the handling of
#' Missing At Random (MAR) data in longitudinal modeling.
#'
#' @details
#' **Clinical Context:**
#' Following a cardiac event, continuous moderate physical activity is crucial
#' to prevent relapse. Patients in this dataset were instructed to track their
#' daily active minutes.
#'
#' **The Missing Data Mechanism (MAR):**
#' Missing data in clinical trials is rarely completely random. In this dataset,
#' dropout is driven by baseline "Cardiac Self-Efficacy" (`efficacy`), which
#' measures a patient's confidence in their heart's resilience. Patients with
#' lower efficacy scores struggle more as the program progresses. Starting around
#' week 5, these patients feel discouraged by their lack of progress and stop
#' transmitting their tracking data (a classic Missing At Random mechanism).
#'
#' If these missing values are ignored (e.g., via Complete Case Analysis), the
#' overall treatment effect will be artificially inflated (Success Bias), as only
#' the highly self-efficacious and active patients remain in the later weeks.
#'
#' **Data Generation (The "Trumpet" Effect):**
#' Patient trajectories exhibit heteroscedasticity (a widening "trumpet" shape
#' over time). While baseline activity is relatively clustered, the individual
#' slopes vary significantly (SD = 2.5), causing the data points to spread out
#' heavily in later weeks right as the missing data mechanism cuts off the lower
#' tail of the distribution.
#'
#' @format A tibble in long format with 1,080 rows and 4 variables:
#' \describe{
#'   \item{id}{Unique patient identification number (1 to 90).}
#'   \item{week}{Measurement time point, representing the week of rehab (1 to 12).}
#'   \item{efficacy}{Baseline Cardiac Self-Efficacy score (0 to 100). Lower scores indicate higher fear-avoidance.}
#'   \item{active_minutes}{Average daily physical activity in minutes. Contains `NA` values for dropped-out observations.}
#' }
#'
#' @source Simulated for pedagogical purposes to demonstrate Missing At Random (MAR) mechanisms.
#'
#' @examples
#' data(cardiac_rehab)
"cardiac_rehab"
