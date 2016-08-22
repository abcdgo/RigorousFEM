#!/bin/bash

NUMBER_OF_MODS=$1

NUMBER_OF_NON_DISSIPATIVE_MODS=$2
CANDIDATE_ORDER=14
PROOF_ORDER=5

HOW_FAR_IN_TIME=37.69911184308 # == 6 * 2 * pi <->  6* period

CANDIDATE_TIME_STEP=0 #if timeStep == 0 then auto time step
PROOF_TIME_STEP=0.00001

NI=$3
NON_AUTONOMOUS_PART=0.1 #sigma
HULL=$4

time=$(date +"%y%j_%T")

SCRIPT_DIR=Scripts
SRC_DIR=Proof
RESULTS_DIR=experiments/$time
BIN_DIR=bin
mkdir -p $RESULTS_DIR
mkdir -p $BIN_DIR

CANDIDATE_FILE=candidate_order_${CANDIDATE_ORDER}_$time.txt
TIME_FILE=times.txt

if [ ! -f $BIN_DIR/find_candidate ];
then
    g++ -O2 $SRC_DIR/candidate.cpp -o $BIN_DIR/find_candidate -I../capd/capdDynSys4/include `../capd_build/bin/capd-config --cflags --libs`
fi

/usr/bin/time -o $RESULTS_DIR/$TIME_FILE -v ./$BIN_DIR/find_candidate $NUMBER_OF_MODS $CANDIDATE_ORDER $NI $NON_AUTONOMOUS_PART $CANDIDATE_TIME_STEP $HOW_FAR_IN_TIME > $RESULTS_DIR/$CANDIDATE_FILE

LINE=`tail -2 $RESULTS_DIR/$CANDIDATE_FILE | head -1`
RTN=`python $SCRIPT_DIR/candidateToRigorous.py $LINE`
IFS='*' read -ra DATA<<< "$RTN"

STARTING_TIME="${DATA[0]}"
INIT_VECTOR="${DATA[1]}"

if [ ! -f $BIN_DIR/inclusion ];
then
    g++ -fpermissive -O2 $SRC_DIR/inclusion.cpp -o $BIN_DIR/inclusion -I../capd/capdDynSys4/include `../capd_build/bin/capd-config --cflags --libs`
fi

if [ ! -f $BIN_DIR/iter ];
then
    g++ -fpermissive -O2 $SRC_DIR/iter2.cpp -o $BIN_DIR/iter -I../capd/capdDynSys4/include `../capd_build/bin/capd-config --cflags --libs`
fi

TMP_FILE1=tmp1.txt
TMP_FILE2=tmp2.txt


    PROOF_FILE=proof_order_${PROOF_ORDER}_$time.txt
    SUMMARY_FILE=summary_order_${PROOF_ORDER}_$time.txt

    /usr/bin/time -a -o $RESULTS_DIR/$TIME_FILE -v ./$BIN_DIR/inclusion $NUMBER_OF_MODS $NUMBER_OF_NON_DISSIPATIVE_MODS $PROOF_ORDER $NI $NON_AUTONOMOUS_PART $PROOF_TIME_STEP $STARTING_TIME $HULL $INIT_VECTOR > $RESULTS_DIR/$PROOF_FILE
    SUBSET=`grep "Subset" $RESULTS_DIR/$PROOF_FILE`

    FFILE=$PROOF_FILE

PERIODS=0
while :
do
    SUBSET=`grep "Subset" $RESULTS_DIR/$FFILE`
    if [ "$SUBSET" = "Subset: 0" ];
    then
        PERIODS=$(( $PERIODS + 1 ))
        INIT_VECTOR_1=`grep -n "Intsec:" $RESULTS_DIR/$FFILE | cut -f1 -d: | tail -1`
        INIT_VECTOR_1=`sed -n "$INIT_VECTOR_1 p" $RESULTS_DIR/$FFILE`
        INIT_VECTOR_1=`echo ${INIT_VECTOR_1#*:}`
        FFILE=proof_iter_${PERIODS}_order_${PROOF_ORDER}_$time.txt
        SUMMARY_FILE=summary_iter_${PERIODS}_order_${PROOF_ORDER}_$time.txt
        /usr/bin/time -a -o $RESULTS_DIR/$TIME_FILE -v ./$BIN_DIR/iter $NUMBER_OF_MODS $NUMBER_OF_NON_DISSIPATIVE_MODS $PROOF_ORDER $NI $NON_AUTONOMOUS_PART $PROOF_TIME_STEP $STARTING_TIME $HULL $INIT_VECTOR_1 > $RESULTS_DIR/$FFILE


        TMP1=`grep -n "Time" $RESULTS_DIR/$FFILE | cut -f1 -d: | head -1`
        TMP2=`grep -n "Time" $RESULTS_DIR/$FFILE | cut -f1 -d: | tail -1`
        sed -n "$TMP1 p" $RESULTS_DIR/$FFILE > $RESULTS_DIR/$TMP_FILE1
        sed -n "$TMP2 p" $RESULTS_DIR/$FFILE >> $RESULTS_DIR/$TMP_FILE1

        echo "Time" > $RESULTS_DIR/$SUMMARY_FILE
        ./$SCRIPT_DIR/checkProof.py $RESULTS_DIR/$TMP_FILE1 >> $RESULTS_DIR/$SUMMARY_FILE

        TMP1=`grep -n "Dbase" $RESULTS_DIR/$FFILE | cut -f1 -d: | head -1`
        TMP2=`grep -n "Dbase" $RESULTS_DIR/$FFILE | cut -f1 -d: | tail -1`
        sed -n "$TMP1 p" $RESULTS_DIR/$FFILE > $RESULTS_DIR/$TMP_FILE2
        sed -n "$TMP2 p" $RESULTS_DIR/$FFILE >> $RESULTS_DIR/$TMP_FILE2

        echo "Dbase" >> $RESULTS_DIR/$SUMMARY_FILE
        ./$SCRIPT_DIR/checkProof.py $RESULTS_DIR/$TMP_FILE2 >> $RESULTS_DIR/$SUMMARY_FILE

        echo "rzad inkluzji: $PROOF_ORDER rzad kandydata: $CANDIDATE_ORDER liczba modow: $NUMBER_OF_MODS liczba niedysypatywnych: $NUMBER_OF_NON_DISSIPATIVE_MODS hull: $HULL ni: $NI non autonomous parameter: $NON_AUTONOMOUS_PART" | cat - $RESULTS_DIR/$SUMMARY_FILE > $RESULTS_DIR/temp && mv $RESULTS_DIR/temp $RESULTS_DIR/$SUMMARY_FILE
    else
        break
    fi

done

    SUMMARY_FILE=summary_order_${PROOF_ORDER}_$time.txt

    TMP1=`grep -n "Time" $RESULTS_DIR/$PROOF_FILE | cut -f1 -d: | head -1`
    TMP2=`grep -n "Time" $RESULTS_DIR/$PROOF_FILE | cut -f1 -d: | tail -1`
    sed -n "$TMP1 p" $RESULTS_DIR/$PROOF_FILE > $RESULTS_DIR/$TMP_FILE1
    sed -n "$TMP2 p" $RESULTS_DIR/$PROOF_FILE >> $RESULTS_DIR/$TMP_FILE1

    echo "Time" > $RESULTS_DIR/$SUMMARY_FILE
    ./$SCRIPT_DIR/checkProof.py $RESULTS_DIR/$TMP_FILE1 >> $RESULTS_DIR/$SUMMARY_FILE

    TMP1=`grep -n "Dbase" $RESULTS_DIR/$PROOF_FILE | cut -f1 -d: | head -1`
    TMP2=`grep -n "Dbase" $RESULTS_DIR/$PROOF_FILE | cut -f1 -d: | tail -1`
    sed -n "$TMP1 p" $RESULTS_DIR/$PROOF_FILE > $RESULTS_DIR/$TMP_FILE2
    sed -n "$TMP2 p" $RESULTS_DIR/$PROOF_FILE >> $RESULTS_DIR/$TMP_FILE2

    echo "Dbase" >> $RESULTS_DIR/$SUMMARY_FILE
    ./$SCRIPT_DIR/checkProof.py $RESULTS_DIR/$TMP_FILE2 >> $RESULTS_DIR/$SUMMARY_FILE

    echo "rzad inkluzji: $PROOF_ORDER rzad kandydata: $CANDIDATE_ORDER liczba modow: $NUMBER_OF_MODS liczba niedysypatywnych: $NUMBER_OF_NON_DISSIPATIVE_MODS hull: $HULL ni: $NI non autonomous parameter: $NON_AUTONOMOUS_PART" | cat - $RESULTS_DIR/$SUMMARY_FILE > $RESULTS_DIR/temp && mv $RESULTS_DIR/temp $RESULTS_DIR/$SUMMARY_FILE

    rm $RESULTS_DIR/*tmp*

