<!--
SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@physik.uni-freiburg.de>
SPDX-FileCopyrightText: 2021 Andreas Härtel <http://andreashaertel.anno1982.de/>

SPDX-License-Identifier: CC-BY-SA-4.0
-->

# git instructions

Tell git that this is a git-repository folder: 
```bash
cd <PROJECTFOLDER>
git init
```

Link the working directory to a remote repository. "origin" is the alias for this remote repository.
```bash
git remote add origin https://github.com/andreashaertel/capdft.git
```

Show the remote repositories linked to this working directory
```bash
git remote -v
```

Tell git to remember your credentials
```bash
git config --global credential.helper store
```

Tell git what your favourite editor is.
```bash
git config --global core.editor vim
git config --global diff.tool vimdiff
```

Get the latest version of the master branch from the remote repository
```bash
git pull origin master
```

Add a file to the staging area (make it ready for commit)
```bash
git add <FILENAME>
```

Revert changes instead of adding them to staging
```bash
git checkout -- <FILENAME>
```

Commit changes in the staging area to your local working directory
```bash
git commit -m "<MESSAGE>"
```

Push changes to a branch (e.g. master) of the remote repository
```bash
git push origin master
```

Make a new branch locally
```bash
git checkout -b <NEWBRANCHNAME>
```

Make a new branch copied from another branch locally
```bash
git checkout -b <NEWBRANCHNAME> <EXISTINGBRANCHNAME>
```

Show all branches
```bash
git branch
```

Merge onto master. First switch to master branch then ...
```bash
git merge <branch>
```

Have a look at this cheat sheet 
[https://rogerdudler.github.io/git-guide/](https://rogerdudler.github.io/git-guide/ "Have a look ...")

Use ssh to save time 
[https://github.com/settings/keys](https://github.com/settings/keys)

