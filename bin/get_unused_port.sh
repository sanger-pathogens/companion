#!/bin/sh

for port in $(seq 4444 65000); do
  echo -ne "\035" | telnet 127.0.0.1 $port > /dev/null 2>&1;
  [ $? -eq 1 ] && echo "$port" && break;
done