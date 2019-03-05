#!/bin/bash

function dotest() {

  lib=$1
  cov=$2
  status=0

  if ! command -v lcov > /dev/null; then
    cov=""
  fi

  echo "Compiling tests against $lib..."
  make clean
  if ! make LIB=$lib > /dev/null; then
    echo "$lib compilation failed!"
    return 1
  fi


  if [ "$cov" != "" ]; then
    lcov -c -i -b . -d . -o cov.baseline > /dev/null

  fi

  echo "Running tests against $lib..."
  for srcfile in *.cpp; do
    execfile=`basename $srcfile .cpp`

    if [ ! -f $execfile ]; then
      continue
    fi

    outfile="$lib"_"$execfile".txt
    benchmark=benchmarks/$outfile

    if ! ./$execfile > $outfile; then
      echo "$execfile[$lib] errored!"
      status=1
      continue
    fi

    if ! cmp -s $outfile $benchmark; then
      echo "$execfile[$lib] failed!"
      status=1
      continue
    fi

    rm $outfile
  done

  if [ "$cov" != "" ]; then
    echo "Generating coverage report..."
    lcov -c -d . -b . -o cov.out                      > /dev/null
    lcov -a cov.baseline -a cov.out -o cov.all        > /dev/null
    lcov --extract cov.all "*armawrap*"   -o cov.all  > /dev/null
    lcov --remove  cov.all "*tests*"      -o cov.all  > /dev/null
    lcov --remove  cov.all "*armadillo*"  -o cov.all  > /dev/null
    report=`lcov --list    cov.all`

    echo "$report"

    ncovered=`echo "$report" | grep ".hpp" | wc -l`
    ntotal=`ls ../armawrap/*.hpp | wc -w`
    echo "$ncovered / $ntotal armawrap files covered"

  fi

  rm -f cov.* *.gcda *.gcno

  return $status
}

pushd `dirname $0` > /dev/null

echo "Compiling newmat ..."
pushd newmat > /dev/null
make clean
make
popd > /dev/null

exitcode=0


if ! dotest newmat ""; then
  exitcode=1
fi


if ! dotest armawrap "1"; then
  exitcode=1
fi

echo "Done!"
make clean
popd > /dev/null

exit $exitcode
