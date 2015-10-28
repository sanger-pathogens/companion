#!/bin/sh

i=1
while [ "$i" -lt 10 ]; do
  port=`awk -v min=4444 -v max=65000 'BEGIN{srand(systime() + PROCINFO["pid"]); print int(min+rand()*(max-min+1))}'`;
  nc -z -w5 127.0.0.1 $port 2>&1;
  [ $? -eq 1 ] && echo "$port" && break;
  >&2 echo "retrying port $port, try $i"
  i=`expr $i + 1`
done