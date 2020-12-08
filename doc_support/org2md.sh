#!/usr/bin/bash

# name:    org2md.sh
# author:  nbehrnd@yahoo.com
# license: MIT, 2020
# date:    2020-12-08 (YYYY-MM-DD)
# edit:
#

# Convert .org documentation into GitHub flavour of markdown
#
# While large sections of this fork are documented in Emacs' org mode,
# the footnote's rendering its markdown export does not fit the
# dialect used by GitHub.  It is easier to let pandoc perform this
# conversion instead.

pandoc -s README.org -t gfm -o README.md

