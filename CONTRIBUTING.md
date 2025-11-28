# Please contribute to pyMBE!

## Bug reports
If you discover any strange or unexpected behaviour in pyMBE,
please report it by opening an issue in this GitHub repository or by contacting
Dr. Pablo M. Blanco ([@pm-blanco](https://github.com/pm-blanco)) or any other member of the pyMBE development team.
Once a ticket is open in GitHub, we will work to fix the issue as soon as possible.

## New contributors
New developers are welcome to contribute to extend the functionalities of the pyMBE module.
To contribute to the pyMBE module, first one needs to be added as a member of this GitHub repository.

If you want to contribute to the development of the pyMBE module,
please contact Dr. Pablo M. Blanco (<pablb@ntnu.no>).

## Authorship and contributorship policy
For more information on our authorship and contributorship policy,
we refer the interested reader to our [Contributor Manifesto](CONTRIBUTOR_MANIFESTO.md).

## How to contribute to pyMBE
Create a fork of the repository and submit your contribution in the form of a pull request.
Give your branch a short and meaningful name, and start your work from the most recent commit on the main branch.
Any new version of the code must reproduce all the data stored in `testsuite/data`.
Before pushing your code, run `make tests` to execute the full test suite,
`make pylint` to check for code issues, and `make docs` to confirm sphinx can build the user guide.
When rapidly prototyping code, run `make unit_tests` periodically to check for regressions.
All new code will be reviewed by at least one member of the pyMBE development team before
merging into the main branch to ensure that a functional version of the code is always available.
Class methods are sorted in alphabetical order.
We follow semantic versioning and keep a changelog.

## Community standards
We aim to keep the pyMBE community an open and friendly space for researchers to collaborate
in developing tools for building coarse-grained models of polyelectrolytes and biomolecules.
We adhere to the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md) v2.1
to provide guidelines on the expected behaviour of the members of our community
and help maintain a healthy and a positive environment for our developers and users.
