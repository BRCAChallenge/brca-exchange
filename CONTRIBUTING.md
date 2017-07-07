# Contributing Guide

If you do not have write permissions, but would like to make a contribution, please refer to the following guidelines.

## Pull Request Guidelines

Before submitting a pull request, please make sure to meet the following:

1. Include tests for bug fixes in server side code. Make sure the tests fail if the bug persists and pass if the bug is resolved.
1. Make any necessary updates to documentation.
1. Resolve all merge conflicts before submitting a PR.


## Suggested Workflow

### First, Create a Fork

1. Visit https://github.com/BRCAChallenge/brca-exchange and click the "Fork" button.
1. Once your fork is created, use your favorite git client to clone the repo to your local machine, or from the command line:

```shell
# Clone your fork to your local machine
git clone git@github.com:USERNAME/brca-exchange.git
```

## Maintain Your Fork

Because your fork is separate from the main repo, you'll want to track the "upstream" repo that you forked. First, add a remote:

```shell
# Add 'upstream' repo to list of remotes
git remote add upstream https://github.com/brcachallenge/brca-exchange.git

# Verify the new remote named 'upstream'
git remote -v
```

Whenever you want to update your fork with the latest upstream changes, you'll need to first fetch the upstream repo's branches and latest commits to bring them into your repository:

```shell
# Fetch from upstream remote
git fetch upstream

# View all branches, including those from upstream
git branch -va
```

Now, checkout your own master branch and merge the upstream repo's master branch:

```shell
# Checkout your master branch and merge upstream
git checkout master
git merge upstream/master
```

If there are no unique commits on the local master branch, git will simply perform a fast-forward. This should always be the case as no commits should be made directly to your master branch. If changes have been made to your local master, make sure to restore local master to match upstream master.

Now, your local master branch is up-to-date with everything modified upstream.

## Doing Your Work

### Create a Branch
Whenever you begin work on a new feature or bugfix, it's important that you create a new branch. Not only is it proper git workflow, but it also keeps your changes organized and separated from the master branch so that you can easily submit and manage multiple pull requests for every task you complete.

To create a new branch and start working on it:

```shell
# Checkout the master branch - you want your new branch to come from master
git checkout master

# Create a new branch named newfeature (give your branch its own simple informative name)
git branch newfeature

# Switch to your new branch
git checkout newfeature
```

Now, go to town hacking away and making whatever changes you want to.

## Submitting a Pull Request

### Cleaning Up Your Work

Prior to submitting your pull request, you might want to do a few things to clean up your branch and make it as simple as possible for the original repo's maintainer to test, accept, and merge your work.

If any commits have been made to the upstream master branch, you should rebase your development branch so that merging it will be a simple fast-forward that won't require any conflict resolution work.

```shell
# Fetch upstream master and merge with your repo's master branch
git fetch upstream
git checkout master
git merge upstream/master

# If there were any new commits, rebase your development branch
git checkout newfeature
git rebase master
```

This is not required, but it may be desirable to squash some of your smaller commits down into a small number of larger more cohesive commits. You can do this with an interactive rebase:

```shell
# Rebase all commits on your development branch
git checkout
git rebase -i master
```

This will open up a text editor where you can specify which commits to squash.

### Submitting

Once you've committed and pushed all of your changes to GitHub, go to the page for your fork on GitHub, select your development branch, and click the pull request button. If you need to make any adjustments to your pull request, just push the updates to GitHub. Your pull request will automatically track the changes on your development branch and update.

### Sources
This guide borrows heavily from the following resources:
* [Chaser324 - Github Forking Gist](https://gist.github.com/Chaser324/ce0505fbed06b947d962)
* [Enzyme Contributing Guide](https://github.com/airbnb/enzyme/blob/master/CONTRIBUTING.md)
