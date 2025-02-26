# Contribution Guidelines

Thank you for considering contributing to the Raven Hydrological Modelling Framework! By participating in this project, you agree to abide by the following terms and conditions.

## Contribution Process

1. Fork the repository and create your branch from `main`.
2. Make your changes.
3. Ensure your code adheres to the project's coding standards.
4. Update documentation if necessary.
5. Submit a pull request.

## License

By contributing to this project, you agree that your contributions will be licensed under the [Artistic License v2.0](https://github.com/CSHS-CWRA/RavenHydroFramework/blob/main/LICENSE).

## Contributor Licensing Agreement (CLA)

Before we can accept your contributions, you need to agree to our Contributor License Agreement (CLA). This agreement helps us ensure that we can freely distribute your contributions as part of the project.

### Summary of the CLA

* Grant of Copyright License: You grant us a perpetual, worldwide, non-exclusive, royalty-free license to reproduce, prepare derivative works of, publicly display, publicly perform, sublicense, and distribute your contributions and derivative works thereof under the Artistic License v2.0.
* Grant of Patent License: You grant us a perpetual, worldwide, non-exclusive, royalty-free patent license to make, have made, use, offer to sell, sell, import, and otherwise transfer your contributions under certain conditions.
* Source of Contribution: You represent that your contributions are your original creation and that you have the right to grant these licenses.
* Rights in Contributions: You confirm that you have the necessary rights and permissions to make your contributions, including any required approvals from your employer.
* Representations and Warranties: You represent that you are legally entitled to grant the above licenses and that your contributions are your original creation.

### Full Agreement

The full text of the CLA can be found [here](https://github.com/CSHS-CWRA/RavenHydroFramework/blob/main/contributor-licensing-agreement.txt).

### How to Agree to the CLA

To agree to the CLA, you must sign the checkbox pertaining to it in your first pull request:

" - [x] I have read and agree to the terms of the Contributor License Agreement."

Alternatively, you can sign the CLA by submitting the completed CLA form to [James Craig](mailto:jrcraig@uwaterloo.ca?subject=RAVEN-CLA).

## Pull Requests

When opening a Pull Request, we ask that contributors label their changes with an appropriate title and provide descriptions of the changes that they wish to implement, whether those changes may cause disruptions to the normal behaviour of the project ("Breaking Changes"), and provide links to any pertinent GitHub Issues, Pull Requests, or other information.

### Styling

This project employs a handful of code styling conventions. Those currently in use are:

* No whitepaces at the end of lines.
* All files must end with a newline ('\n').

Coding standards are enforced via [pre-commit](https://pre-commit.com/), with automatic corrections performed via [pre-commit CI](https://pre-commit.ci/).

### Documentation

Documentation for the Raven Hydrological Modelling Framework can be found at https://raven.uwaterloo.ca/Downloads.html.

### Testing

This library does not currently use testing. It is expected that all new contributions should compile for the following architectures:

- Linux (x86_64)
- Windows (x86_64)
- macOS (x86_64, ARM64)
