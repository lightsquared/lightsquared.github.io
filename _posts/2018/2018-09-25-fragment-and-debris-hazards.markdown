---
title: Fragment-and-Debris-Hazards
date: '2018-09-25 11:03'
categories: engineering
tags:
  - hazard
  - safety
  - fragment
  - debris
---

When the [United States Department of Defense Explosives Safety Board (DDESB)][9615d5d6] determines fragment and debris hazards they use a 6-step process based on [Technical Paper 12][1096b4d4] (TP-12).  I was using this process recently in a paper and ran down a rabbit hole when I got to step 4.  Before I drag you down the hole, let me paraphrase *all* the steps as they are listed in TP-12.

1.  1. Determine initial fragment velocity either using the appropriate Gurney equations (or determine experimentally). The general Gurney equation is given by:

$$v_0= \frac{\sqrt{2E}}{\frac{M}{C}+\frac{n}{n+2}}$$

where $\sqrt{2E}$ is the Gurney energy constant, $M$ is the mass of the explosive, $C$ is the mass of the casing, and $n=1,\: 2,$ or $3$ for plane, cylindrical, and spherical symmetry.

2.  Estimate the average fragment mass $m_0$ by fitting the Mott distribution to a single-weapon arena test.  Where the Mott distribution is given by:

$$N=\frac{M_T}{m_0}e^{-\left(\frac{2m}{m_0}\right)^{\frac{1}{2}}}$$

where $N$ is the number of fragments heavier than mass $m$ and $M_T$ is the total mass of all the fragments.

3.  Determine the mass $m$ of the lightest hazardous fragment reaching a specified distance $R$ either from the solution of:

$$2E_{cr}=mV^2e^{\left( \frac{-2R}{L_1m^{\frac{1}{3}}}\right)}$$
or from the solution of:

$$2E_{cr}=gL_1m^{\frac{4}{3}}$$
whichever gives the smaller value of $m$.  The critical level of kinetic energy at impact is defined as $58\:ft-lb\: \left( 79\:J \right)$ and $L_1$ is defined as:
$$L_1=\frac{2\left(k^2m\right)^{\frac{1}{2}}}{C_D\rho}$$

where $m$ is the fragment mass, $k$ is a shape factor, and $C_D$ is the coefficient of drag.

5.  Determine the areal density of fragments heavier than $m$ reaching distance $R$ from the inverse-square law:

$$$$

  [1096b4d4]: http://www.esd.whs.mil/Portals/54/Documents/FOID/Reading%20Room/Other/10-F-0806_Fragment_and_Debris_Hazards.pdf "Fragment and Debris Hazards"
  [9615d5d6]: https://www.denix.osd.mil/ddes/home/ "DDESB"
