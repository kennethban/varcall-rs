#!/usr/bin/env bash
# A truth case's arguments (truth/<case>.cmd). Paths into the data directory are stored as $VARCALL_DATA/..., so the
# truth folders work wherever the data directory is.
#   cmdargs.sh FILE        FILE's arguments with $VARCALL_DATA expanded (readers: $($HARNESS/cmdargs.sh truth/$c.cmd))
#   echo ARGS | cmdargs.sh -w > FILE   ARGS with the data directory's path written as $VARCALL_DATA (writers)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
d=${VARCALL_DATA%/}
if [ "$1" = -w ]; then
  p=$(cd "$d" && pwd -P)
  sed -e "s|$p/|\$VARCALL_DATA/|g" -e "s|$d/|\$VARCALL_DATA/|g"
else
  sed "s|[\$]VARCALL_DATA|$d|g" "$1"
fi
