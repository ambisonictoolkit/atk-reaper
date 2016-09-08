#!/bin/sh

if [ -d ~/Library/Application\ Support/REAPER/Effects/ATK ]; then
  rm -rf ~/Library/Application\ Support/REAPER/Effects/ATK
fi

exit 0


#[ -d ~/Library/Application\ Support/REAPER/Effects/ATK ] && echo 'Directory found' || echo 'Directory /tmp not found'

