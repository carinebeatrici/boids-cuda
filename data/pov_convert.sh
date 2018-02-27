#!/bin/bash
for i in *.pov
do
  povray $i +W800 +H600 -D 2>/dev/null;
  echo $i
done
