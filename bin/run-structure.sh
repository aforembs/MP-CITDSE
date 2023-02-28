#!/bin/bash

function usage {
        echo "Usage: $(basename $0) [-he:f:i:]" 2>&1
        echo '   -h   shows this help message'
        echo '   -e   the numberof electrons 1 or 2'
        echo '   -f   path to the .yaml input file'
        echo '   -i   path to directory of the cfg-<L>.inp files only used if -e 2'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring=":he:f:i:"

while getopts ${optstring} arg; do
  case "${arg}" in
    h)  usage 
        ;;
    e) NUM_E=${OPTARG} ;;
    f) YAML_FILE=${OPTARG} ;;
    i) CFG_DIR=${OPTARG} ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

if [[ ${NUM_E} == 1 ]]; then
  bin/h1e -f ${YAML_FILE}
  bin/w1e -f ${YAML_FILE}
  bin/d1e -f ${YAML_FILE}
elif [[ ${NUM_E} == 2 ]]; then
  bin/h1e -f ${YAML_FILE}
  bin/w1e -f ${YAML_FILE}
  bin/d1e -f ${YAML_FILE}
  bin/id2ec -f ${YAML_FILE} -i ${CFG_DIR}
  bin/r12 -f ${YAML_FILE} -i ${CFG_DIR}
  bin/d2e -f ${YAML_FILE} -i ${CFG_DIR}
  bin/w2e -f ${YAML_FILE}
else
  echo "Invalid number of electrons: ${NUM_E}."
  echo "Only 1 or 2 supported."
fi