#! /bin/sh

# generate ChangeLog from git shortlog
echo "rebuilding ChangeLog"
# loop through sequential tags to build ChangeLog
# 0.0.1..0.0.2
# 0.0.2..0.0.3
# etc
[[ -f ChangeLog ]] && rm ChangeLog
T=( $(git tag -l) )
for i in $(seq 1 $((${#T[@]}-1)) | tac); do
  echo -ne "Release ${T[$i]}:                                          " >> ChangeLog
  date --date=- +'%b %d, %Y' --date="$(git show -s --format="%ci" ${T[$i]} | tail -1)" >> ChangeLog
  git shortlog -w79,6,9 ${T[$i-1]}..${T[$i]} | grep -v 'minor:' | grep -e '^ ' >> ChangeLog
  echo "" >> ChangeLog
done

if [[ ! -f COPYING ]]; then
  echo "get GPLv3"
  wget -q http://www.gnu.org/licenses/gpl-3.0.txt -O COPYING
fi

if [[ ! -f COPYING.LESSER ]]; then
  echo "get LGPLv3"
  wget -q http://www.gnu.org/licenses/lgpl-3.0.txt -O COPYING.LESSER
fi

echo "update configure.ac, etc"
autoreconf --install --warnings=all $*
