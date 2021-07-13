# Contributing to varCA
Thank you so much for taking the time to contribute to varCA! :rocket::tada:

Contributions are always welcome, and they are greatly appreciated!

## Types of contributions
### Report a bug
If you have found a bug, please report it on [our issues page](issues) rather than emailing us directly. Others may have the same issue and this helps us get that information to them.

Before you submit a bug, please search through our issues to ensure it hasn't already been reported.

The most helpful Github issues include
- the version of varCA you are using, although it's best to use the latest version
- the version of Snakemake you are using
- detailed steps to help us reproduce your error, ideally with the example dataset distributed with varCA
### Fix a bug
Look through our issues page for bugs. We especially need help with bugs labeled "help wanted". If you want to start working on a bug then please write a message within the thread for that issue on our issues page, so that no one is duplicating work.
### Implement a new feature
Our issues page will almost always have features on our wishlist. Once again, if you want to start working on a feature then please write a message within the thread for that feature on our issues page, so that no one is duplicating work.

Have an idea for a new feature that isn't on our wishlist? We'd love to hear about it! Before getting to work, please create a Github issue for it, so that you can make sure we're in agreement about what it should do.

## How to fix a bug or implement a new feature
Please create a pull request! A PR is a collection of changes that you have made to the code that we can review and potentially integrate into varCA.

To create a pull request you need to do these steps:

1. Create a Github account.
2. [Fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo#forking-a-repository) the repository.
    - Click the "Fork" button in the top, right corner
    - Or, if you had already forked the repository a while ago, [sync your fork](https://docs.github.com/en/github/collaborating-with-pull-requests/working-with-forks/syncing-a-fork) to make sure you're working with the latest version of varCA.
4. [Clone your fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo#cloning-your-forked-repository) locally.
5. `cd varCA` into the new directory
6. Create a new branch with `git checkout -b <descriptive_branch_name>`
7. Make your changes to the code.
8. Test that they work. And test your code with any existing tests in the repository to ensure you haven't broken anything. See "Testing" below.
9. Please add any comments to the documentation that would help users understand how to use your new code.
10. Commit your changes. Please use informative commit messages and do your best to ensure the commit history is clean and easy to interpret.
11. Now you can push your changes to your Github copy of varCA by running `git push origin <descriptive_branch_name>`
12. Go to your Github copy of varCA in your browser and create a pull request. Be sure to change the pull request target branch to `master` on this original repository!
13. Please write an informative pull request detailing the changes you have made and why you made them. Tag any related issues by referring to them by a hashtag followed by their ID.

### Testing
It's critical that you test your new code to make sure that you're not inadvertantly introducing any new bugs.
If you are fixing a bug or implementing a new feature, please add at least one test to cover your new code.
In addition, please also test the entire pipeline with the example dataset to ensure existing code hasn't broken.

## Style
### Code
- Please use 4 spaces instead of tabs
- Please adhere to PEP8 whenever possible
- Try not to exceed 88 characters per line
- Do not duplicate strings within Snakefiles. Please refer to them as `rules.rulename.output.outputname` whenever possible.
### Git Commit Messages
- Use the present tense ("Add feature" not "Added feature")
- Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
- Reference issues and pull requests liberally after the first line
