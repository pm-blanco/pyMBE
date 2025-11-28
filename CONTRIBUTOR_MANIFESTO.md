# Manifesto for Open Source Software Authorship and Credit

This document establishes our vision of open source software contributorship and credit in scholarly publishing.
We follow the structure outlined in section "Establish an authorship policy" of Garijo et al. 2022[^Garijo2022]
and re-use core elements of "SORTÆD: Software Role Taxonomy and Authorship Definition" version 0.1[^Leem2023].

We welcome and value all contributions, from small fixes to large features.
We also recognize the importance of non-code contributions, including, but not limited to:
bug reports, code reviews, admin tasks, wiki contributions, and feedback from users and scientific advisors.

We agree with the CRediT's assessment that traditional bibliographic conventions for defining authorship are outdated[^ANSI_NISO_2022].
We agree that a new form of credit needs to emerge[^Brand2015] that would allow us to move from authorship to contributorship.
We also think the current CRediT taxonomy is inadequate to describe software contributions.

We believe that adopting a generic name such as "pyMBE Team", "Project Jupyter",
or "TensorFlow Developers" tends to invisibilize individual contributors.
We believe that listing all contributors as co-authors regardless
of their level of contribution dilutes the value of authorship,
and may be detrimental to authors who have long left the software community
and no longer need (or desire) to be associated with it in scholarly publishing.
We believe credit should go to contributors who significantly influenced the software.
We understand and accept the inherent subjectivity of determining what constitutes a significant contribution.
We believe this assessment needs to be revised for every major release of the software.
Major releases are defined in accordance with Semantic Versioning 2.0.0[^SemVer].

We make a distinction between authorship w.r.t. intellectual property,
for which any code contribution grants authorship as per copyright law
(and when applicable per a Contributor License Agreement),
versus authorship w.r.t. the software project,
for which code contributions aren't necessary to qualify.
For example, a code reviewer who has a significant impact on the software
can be credited as an author of the software in scholarly publishing,
but cannot claim ownership of the source code.
Inversely, a code contributor who made trivial changes to the source code
that din't impact the software in a meaningful way won't necessarily
be eligible for credit in scholarly publishing, but can still claim
ownership as defined in copyright law.[^Pamela2018]

In the following, we'll use the term "contributor" to reflect any individual
who made code contributions and non-code contributions to the software.
This definition is completely orthogonal to the definition of "contributor"
given in the license agreement, which is mainly concerned with intellectual
property and redistribution rights.
The present document is only concerned with scholarly publishing,
and how credit can be given to the software contributors.
We specifically use this term to detach it from the traditional "author" term,
which in our view has a narrower definition which naturally excludes non-code contributors.
We will use the term "author" to refer to the author list in a citation.
Our users are expected to acknowledge the use of the software with at least two citations[^Smith2016]:
first is the software paper, whose author list was determined according to the CRediT taxonomy,
second is the software source code, which can be either a stable release of the
software with semantic versioning, or a development version with a commit hash;
in both cases, the author list of the software corresponds to the contributor
list as defined in this document, and is encoded in the citation file format.[^Druskat2021]

We believe that significant contributions need to be credited. We also recognize
that historic contributors who are no longer active in the project by choice
may no longer be interested in being prominently credited for the software
in scholarly publishing, and cannot be expected to vouch for the current state
and future directions of the project.
In addition, new contributors might benefit from being in the spotlight for their career
development, while established researchers rely on other metrics to advance their career.
Assessing credit on a per-release basis allows us to keep the author list short when
citing the software, while rewarding new contributors for their investment in the project.
This workflow is not unlike the process of determining the author list
on a software paper using the CRediT taxonomy.
When a list of contributors is determined for a software release X.Y.0,
all future bugfix releases X.Y.1 to X.Y.n should feature the same contributors,
just like a software paper errata should feature the same authors.

What constitutes a significant contribution is necessarily subjective,
and cannot be adequately captured by key performance indicators.
We think that developing new features or extending existing features
are excellent ways of securing credit in the next release of the software.
Maintenance efforts should also be rewarded, since a software that can no longer
be installed on a modern operating system due to broken dependencies will quickly fade into irrelevance.

Following these considerations, we recognize the following as valid contributions warranting the status of _contributor_:

- Extension of the functionalities of the software.
- Bugfixes and enhancements of the code quality, following best coding practices.
- Development of tutorials, sample scripts and other documentation of the software.
- Ideas on software design and architecture. This includes ideas on improving the code quality,
  existing features of the software, adding new features and other conceptual contributions
  that have a significant impact on the development of the software.

We classify our contributors as either _active_ contributors or _past_ contributors:

- Active contributors are those who have contributed to either the current
  stable version of the software or the one in development on the main branch.
   - Active contributors are responsible for maintaining the software and
     addressing bug reports on the main branch and last stable release.
   - Only active contributors are eligible for co-authorship on stable releases
     of the software and in the `CITATION.cff` file of the main branch.
   - Active contributors may be invited in prospective scholarly publications,
     depending on their availability and interest to participate in the preparation of the publication.
- Past contributors are those who significantly contributed to any past stable version of the software
  but have not contributed to the code in the last stable version or in the main branch.
   - Past contributors are not responsible for the maintenance of the software.
   - Past contributors are acknowledged as authors of the software for legal purposes,
     but they are not eligible for co-authorship on stable releases
     of the software and in the `CITATION.cff` file of the main branch.

We distinguish the following different roles for active contributors:

- Admin contributors are mainly responsible for the maintenance of the software,
  management of the software community and integration of new features to the software.
    - Admins must participate in the review process of pull requests.
    - Admins will be listed as co-authors on stable releases of the software
      and in the `CITATION.cff` file.
    - The admin team has the final word in decisions concerning major developments of the software.
- Tier-0 contributors are developers of key new features of the current version of the software
  or whose ideas had a substantial impact on the development of the current version of the software.
    - Tier-0 contributors can be occasionally requested to review pull requests.
    - Tier-0 contributors are responsible to develop proper unit and functional tests for new features.
    - Tier-0 contributors will be listed as co-authors on stable releases of the
      software and in the `CITATION.cff` file.
- Tier-1 contributors are occasional developers of individual methods of the software,
  or developers who contributed with  bugfixes and other small code enhancements.
    - Tier-1 contributors may be occasionally invited to review pull requests
      but are under no obligation to accept an invitation.
    - Tier-1 contributors are not obliged to provide functional tests for the
      software but may be requested to contribute to unit tests.
    - Tier-1 contributors will not be listed as co-authors on stable releases
      of the software and in the `CITATION.cff` file.

These roles are not static and will be re-evaluated by the Admin team before every major release of the software.
This review will rely on version control history and the proceedings of our periodic community meetings.
Based on this review, the admin team will propose new roles for active contributors.
These proposals will be communicated to the contributors.
The following scenarios are contemplated:

- If the contributor agrees with the proposed role, the change will take effect in the current major release.
- If the contributor disagrees, the change will not be applied in the current
  major release, but will be re-evaluated in the next major release.
  If the same role is identified during the re-evaluation,
  the new role will be effective in the next major release.
- If no response is received, the Admin team will assume the contributor agrees,
  and the change will be implemented in the next major release.
- Upon request, the Admin team may exceptionally assess and update the role
  of a contributor before a major release in the development version of the software.

In case of conflict between two members of the community, one or more members
of the Admin team not involved in the conflict will mediate and negotiate
with all parties involved to search for a consensus solution.
If all Admin members are deemed non-neutral in the conflict,
then the Admin team will search for an external mediator.
If no consensus is reached after mediation, the Admin team may decide
what they deem best for the development of the software and for a healthy
community environment.

**Bibliography**:

[^Leem2023]: Leem D., Turon G., Gruson H., Chue Hong N., Kaur Bhogal S., Lo S., Druskat S., Soiland-Reyes S., "SORTÆD: Software Role Taxonomy and Authorship Definition," version 0.1. *Zenodo* (2023). doi:[10.5281/zenodo.7896456](https://doi.org/10.5281/zenodo.7896456) and [live version](https://sdruskat.net/software-authorship)
[^ANSI_NISO_2022]: ANSI, NISO, "CRediT, Contributor Roles Taxonomy." ANSI/NISO Standard Z39.104-2022, Feb. 2022. doi:[10.3789/ansi.niso.z39.104-2022](https://doi.org/10.3789/ansi.niso.z39.104-2022)
[^Brand2015]: Brand A., Allen L., Altman M., Hlava M. and Scott J., "Beyond authorship: attribution, contribution, collaboration, and credit." *Learned Publishing*, 28:151-155 (2015). doi:[10.1087/20150211](https://doi.org/10.1087/20150211)
[^SemVer]: Semantic Versioning 2.0.0, <https://semver.org>
[^Garijo2022]: Garijo D, Ménager H, Hwang L, Trisovic A, Hucka M, Morrell T, Allen A, "Task Force on Best Practices for Software Registries, SciCodes Consortium. Nine best practices for research software registries and repositories." *PeerJ Computer Science* 8:e1023 (2022). doi:[10.7717/peerj-cs.1023](https://doi.org/10.7717/peerj-cs.1023)
[^Pamela2018]: Pamela S. Chestek, "A Theory of Joint Authorship for Free and Open Source Software Projects." *Colorado Technology Law Journal* 16(2):285 (2018). <https://scholar.law.colorado.edu/ctlj/vol16/iss2/5/>
[^Smith2016]: Smith A. M., Katz D. S., Niemeyer K. E., FORCE11 Software Citation Working Group, "Software citation principles." *PeerJ Computer Science* 2:e86 (2016). doi:[10.7717/peerj-cs.86](https://doi.org/10.7717/peerj-cs.86)
[^Druskat2021]: Druskat S., Spaaks J. H., Chue Hong N., Haines R., Baker J., Bliven S., Willighagen E., Pérez-Suárez D., Konovalov A., "Citation File Format (1.2.0)." *Zenodo* (2021). doi:[10.5281/zenodo.5171937](https://doi.org/10.5281/zenodo.5171937)
