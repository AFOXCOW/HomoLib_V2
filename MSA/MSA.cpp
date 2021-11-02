#include "MSA.h"
using namespace std;

//the Input_Similar_Table must match the Seq_Element exactly
//for example,   A  T
//            A  2 -1
//            T -1  3
//you can input Input_Similar_Table=[[2,-1],[-1,3]] and the Seq_Element must be "AT".
//you can also input Input_Similar_Table=[[3,-1],[-1,2]] and the Seq_Element must be "TA".
MSA::MSA(vector<string> Input_Sequence, vector<vector<double>> Input_Similar_Table,string Seq_Element) {
    //eg: For protein sequence, the element number is 20 for 20 different amino acid
    Number_of_Seq_Element_Plus_Spaceholder = Seq_Element.size();

    number_of_Sequences = Input_Sequence.size();

    // The number of different aligned group.It will change during iteration aligning.
    // at first,each sequence stand for a group.
    Number_of_not_aligned_group = number_of_Sequences;

	for (int i = 0; i < Input_Sequence.size(); ++i) {
		Sequence_labeled temp;
		temp.sequence = Input_Sequence[i];
		temp.label = i;
        if (temp.sequence.size() > Sequence_MaxLength) {
            Sequence_MaxLength = temp.sequence.size();
        }
		Aligning_Sequence.push_back(temp);
	}
    Sequence_MaxLength += 2000;

    // store the Input_Similar_Table as object`s variable.
	for ( int i = 0; i < Number_of_Seq_Element_Plus_Spaceholder; ++i) { //20 amino acid and 1 spaceholder
		for ( int j = 0; j < Number_of_Seq_Element_Plus_Spaceholder; ++j) {
			Similar_Table[i][j] = Input_Similar_Table[i][j];
		}
	}

    // initialize alphaBetHashTable
    // it store the projection between sequence element and its order in sequence element.
    for (int i = 0; i < Seq_Element.size(); ++i) {
        alphaBetHashTable[int(Seq_Element[i])] = i;
    }
    
}
MSA::~MSA(void) {

}
/*
void MSA::Update_Number_of_not_aligned_group(){
    int temp = 0;
    set<int>labels;
    for(int i=0; i < number_of_Sequences; ++i) {
        labels.insert(Aligning_Sequence[i].label);
    }
    Number_of_not_aligned_group = labels.size();
}
*/

//todo 
//what if the label greater than Number_of_not_aligned_group
void MSA::ChooseAlignedPair() {
    double alignScore = -10e5;
    double midScore = 0;
    for (int i = 0; i < Number_of_not_aligned_group; ++i) {
        for (int j = i + 1; j < Number_of_not_aligned_group; ++j)
        {
            midScore = ScoreAlignAlignments(i, j);
            if (midScore > alignScore)
            {
                alignScore = midScore;
                most_similar_group_one = i;
                most_similar_group_two = j;          
            }
        }
    }
}


int MSA::AlignAndUpdate(int firstSet, int secondSet)
{
    int currentSetOneLength = 0;
    int currentSetTwoLength = 0;
    for (int i = 0; i < number_of_Sequences; ++i) {
        if (Aligning_Sequence[i].label == firstSet) {
            currentSetOneLength = Aligning_Sequence[i].sequence.size();
            break;
        }
    }
    for (int i = 0; i < number_of_Sequences; ++i) {
        if (Aligning_Sequence[i].label == secondSet) {
            currentSetTwoLength = Aligning_Sequence[i].sequence.size();
            break;
        }
    }
    vector<vector<double> > scoreMatrix;
    scoreMatrix.resize(currentSetOneLength + 1);
    for (int i = 0; i < currentSetOneLength + 1; ++i) {
        scoreMatrix[i].resize(currentSetTwoLength + 1);
    }
    scoreMatrix[0][0] = 0;
    for (int i = 1; i <= currentSetOneLength; ++i)
        scoreMatrix[i][0] = scoreMatrix[i - 1][0] + AlignColumnToSpace(firstSet, i);
    for (int j = 1; j <= currentSetTwoLength; ++j)
        scoreMatrix[0][j] = scoreMatrix[0][j - 1] + AlignColumnToSpace(secondSet, j);
    double matchScore, delA, insB;


    for (int i = 1; i <= currentSetOneLength; ++i)
        for (int j = 1; j <= currentSetTwoLength; ++j)
        {
            matchScore = scoreMatrix[i - 1][j - 1] + AlignColumnToColumn(firstSet, i, secondSet, j);
            delA = scoreMatrix[i - 1][j] + AlignColumnToSpace(firstSet, i);
            insB = scoreMatrix[i][j - 1] + AlignColumnToSpace(secondSet, j);
            scoreMatrix[i][j] = MaxofThree(matchScore, delA, insB);
        }

    //We traceback the alignment
    vector<int>tracebackLabel(Sequence_MaxLength, 0xff);
    int alignmentSize = 0;
    //match 0; deletion of A, 1; insertion of B,2
    int i = currentSetOneLength;
    int j = currentSetTwoLength;
    while ((i != 0) || (j != 0)){
        if (i == 0){
            for (j = j; j > 0; j--){
                tracebackLabel[alignmentSize] = 2;
                alignmentSize++;
            }
        }
        if (j == 0){
            for (i = i; i > 0; i--){
                tracebackLabel[alignmentSize] = 1;
                alignmentSize++;
            }
        }
        if ((i != 0) && (j != 0)){
            matchScore = scoreMatrix[i - 1][j - 1] + AlignColumnToColumn(firstSet, i, secondSet, j);
            if (fabs(matchScore - scoreMatrix[i][j]) < 1e-10){
                tracebackLabel[alignmentSize] = 0;
                alignmentSize++;
                i--;
                j--;
                continue;
            }
            delA = scoreMatrix[i - 1][j] + AlignColumnToSpace(firstSet, i);
            if (fabs(delA - scoreMatrix[i][j]) < 1e-10){
                tracebackLabel[alignmentSize] = 1;
                alignmentSize++;
                i--;
                continue;
            }
            insB = scoreMatrix[i][j - 1] + AlignColumnToSpace(secondSet, j);
            if (fabs(insB - scoreMatrix[i][j]) < 1e-10){
                tracebackLabel[alignmentSize] = 2;
                alignmentSize++;
                j--;
                continue;
            }

        }
    }
    string alignmentCopy;
    int charLabel;
    //After Trace, we need to update all the relevant things.
    for (i = 0; i < number_of_Sequences; ++i)
    {
        //update alignment pattern for group one
        if (Aligning_Sequence[i].label == firstSet)
        {
            charLabel = 0;
            for (j = alignmentSize - 1; j >= 0; j--)
            {
                if (tracebackLabel[j] == 0 || tracebackLabel[j] == 1)
                {
                    alignmentCopy.push_back(Aligning_Sequence[i].sequence[charLabel]);
                    charLabel++;
                }
                else
                {
                    alignmentCopy.push_back('-');
                }
            }
            Aligning_Sequence[i].sequence = alignmentCopy;
        }
        alignmentCopy.clear();

        //update alignment pattern for group two
        if (Aligning_Sequence[i].label == secondSet)
        {
            charLabel = 0;
            for (j = alignmentSize - 1; j >= 0; j--)
            {
                if (tracebackLabel[j] == 0 || tracebackLabel[j] == 2)
                {
                    alignmentCopy.push_back(Aligning_Sequence[i].sequence[charLabel]);
                    charLabel++;
                }
                else
                {
                    alignmentCopy.push_back('-');
                }
            }
            Aligning_Sequence[i].sequence = alignmentCopy;
        }
        alignmentCopy.clear();
    }
    // new group label the last.
    // eg: old 0~9, after alignment, new group label 8,and others 0~7 whatever.
    // update the sequences` label.

    //record the old group and its label.
    vector<int>interMidSetLabel(sequenceMaxNumber, 0xff);

    //distinguish the new "Number_of_not_aligned_group - 2" and old "Number_of_not_aligned_group - 2"
    //0 for new label, others for old label.
    vector<int>currentLabel(sequenceMaxNumber, 0xff); 
    for (i = 0; i < number_of_Sequences; ++i)
    {
        if (Aligning_Sequence[i].label == firstSet || Aligning_Sequence[i].label == secondSet)
        {
            Aligning_Sequence[i].label = Number_of_not_aligned_group - 2; //new group label.
            currentLabel[i] = 0;
        }
        else
        {   //record the old label and value=0 for old group.
            interMidSetLabel[Aligning_Sequence[i].label] = 0;
        }
    }
    int midValue = 0;
    int midUpdate = 0;//new label start with 0
    for (i = 0; i < Number_of_not_aligned_group; ++i)
    {
        if (interMidSetLabel[i] == 0)//old group
        {
            midValue = i;//old label

            for (j = 0; j < number_of_Sequences; ++j)
            {
                if (Aligning_Sequence[j].label == midValue && currentLabel[j] != 0)//old group with same label.
                {
                    Aligning_Sequence[j].label = midUpdate;
                    currentLabel[j] = 0;
                }
            }
            midUpdate++;
        }
    }
    Number_of_not_aligned_group = Number_of_not_aligned_group - 1;
    scoreMatrix.clear();
    return 0;
}

//two sequences alignment, and score it.
//something like Needleman-Wunsch alg.
double MSA::ScoreAlignAlignments(unsigned int firstSet, unsigned int secondSet)
{
    int currentSetOneLength=0;
    int currentSetTwoLength=0;
    for (int i = 0; i < number_of_Sequences; ++i) {
        if (Aligning_Sequence[i].label == firstSet){
            currentSetOneLength = Aligning_Sequence[i].sequence.size();
            break;
        }
    }
    for (int i = 0; i < number_of_Sequences; ++i) {
        if (Aligning_Sequence[i].label == secondSet){
            currentSetTwoLength = Aligning_Sequence[i].sequence.size();
            break;
        }
    }

    //initialize score Matrix
    vector<vector<double> > scoreMatrix;
    scoreMatrix.resize(currentSetOneLength + 1);
    for (int i = 0; i < currentSetOneLength + 1; ++i) {
        scoreMatrix[i].resize(currentSetTwoLength + 1);
    }
    scoreMatrix[0][0] = 0;


    for (int i = 1; i <= currentSetOneLength; ++i) {
        scoreMatrix[i][0] = scoreMatrix[i - 1][0] + AlignColumnToSpace(firstSet, i);
    }
    for (int j = 1; j <= currentSetTwoLength; ++j) {
        scoreMatrix[0][j] = scoreMatrix[0][j - 1] + AlignColumnToSpace(secondSet, j);
    }
    //A is the vertical axis, B is the horizontal axis
    double matchScore = 0;
    double delA = 0;
    double insB = 0;

    for (int i = 1; i <= currentSetOneLength; ++i) {
        for (int j = 1; j <= currentSetTwoLength; ++j)
        {
            matchScore = scoreMatrix[i - 1][j - 1] + AlignColumnToColumn(firstSet, i, secondSet, j);
            delA = scoreMatrix[i - 1][j] + AlignColumnToSpace(firstSet, i);
            insB = scoreMatrix[i][j - 1] + AlignColumnToSpace(secondSet, j);

            scoreMatrix[i][j] = MaxofThree(matchScore, delA, insB);
        }
    }
    double temp = scoreMatrix[currentSetOneLength][currentSetTwoLength];
    scoreMatrix.clear();
    return temp;
}

// For Group 'setLabel', At position 'columnLocation',when match to the space '-'
// calculate the score for this group.
double MSA::AlignColumnToSpace(int setLabel, int columnLocation)
{   
    int setSequenceTotal = 0;
    char midC = '0';
    int midLabel = 0;
    double returnScore = 0;
    int spaceLabel = 0;

    for (int i = 0; i < number_of_Sequences; ++i) {
        if (Aligning_Sequence[i].label == setLabel) {
            setSequenceTotal++;
        }
    }
    
    vector<double>columnProfile(Number_of_Seq_Element_Plus_Spaceholder,0);

    for (int i = 0; i < number_of_Sequences; ++i) {
        if (Aligning_Sequence[i].label == setLabel){
            midC = Aligning_Sequence[i].sequence[columnLocation - 1];
            midLabel = alphaBetHashTable[int(midC)];
            columnProfile[midLabel]++;
        }
    }
    for (int i = 0; i < Number_of_Seq_Element_Plus_Spaceholder; ++i) {
        columnProfile[i] = columnProfile[i] / setSequenceTotal;
    }

    spaceLabel = alphaBetHashTable[int('-')];
    for (int i = 0; i < Number_of_Seq_Element_Plus_Spaceholder; ++i) {
        returnScore = returnScore + Similar_Table[i][spaceLabel] * columnProfile[i];
    }
    return returnScore;
}

// for group 'setlabelver', at position 'verColumnLocation' 
// for group 'setLabelHor', at position 'horColumnLocation' 
// calculate the score for group.
double MSA::AlignColumnToColumn(int setLabelVer, int verColumnLocation, int setLabelHor, int horColumnLocation)
{
    int setVerSequenceTotal = 0;
    int setHorSequenceTotal = 0;
    char midC = '0';
    int midLabel = 0;
    for (int i = 0; i < number_of_Sequences; ++i) {
        if (Aligning_Sequence[i].label == setLabelVer) {
            setVerSequenceTotal++;
        }
    }
    for (int j = 0; j < number_of_Sequences; ++j) {
        if (Aligning_Sequence[j].label == setLabelHor) {
            setHorSequenceTotal++;
        }
    }
    vector<double> verColumnProfile(Number_of_Seq_Element_Plus_Spaceholder,0);
    vector<double> horColumnProfile(Number_of_Seq_Element_Plus_Spaceholder,0);
    for (int i = 0; i < number_of_Sequences; ++i) {
        if (Aligning_Sequence[i].label == setLabelVer){
            midC = Aligning_Sequence[i].sequence[verColumnLocation - 1];
            midLabel = alphaBetHashTable[int(midC)];
            verColumnProfile[midLabel]++;
        }
    }
    for (int i = 0; i < Number_of_Seq_Element_Plus_Spaceholder; ++i) {
        verColumnProfile[i] = verColumnProfile[i] / setVerSequenceTotal;
    }

    for (int i = 0; i < number_of_Sequences; ++i) {
        if (Aligning_Sequence[i].label == setLabelHor) {
            midC = Aligning_Sequence[i].sequence[horColumnLocation - 1];
            midLabel = alphaBetHashTable[int(midC)];
            horColumnProfile[midLabel]++;
        }
    }

    for (int i = 0; i < Number_of_Seq_Element_Plus_Spaceholder; ++i) {
        horColumnProfile[i] = horColumnProfile[i] / setHorSequenceTotal;
    }
    
    double returnScore = 0;
    for (int i = 0; i < Number_of_Seq_Element_Plus_Spaceholder; ++i) {
        for (int j = 0; j < Number_of_Seq_Element_Plus_Spaceholder; ++j) {
            returnScore = returnScore + Similar_Table[i][j] * verColumnProfile[i] * horColumnProfile[j];
        }
    }
    return returnScore;
}

double MSA::MaxofThree(double A, double B, double C)
{
    double returnValue = A;
    if (B > returnValue) {
        returnValue = B;
    }
    if (C > returnValue) {
        returnValue = C;
    }
    return returnValue;
}

int MSA::IterationAlignProcess()
{
    most_similar_group_one = -1;
    most_similar_group_two = -1;
    for (int i = 1; i <= number_of_Sequences - 1; ++i)
    {
        printf("\n\nIn the %d-th Iteration:\n", i);
        ChooseAlignedPair();
        printf("We Will Align:\n");
        for (int j = 0; j < number_of_Sequences; ++j) {
            if (most_similar_group_one == Aligning_Sequence[j].label) printf("Sequence-%d  ", j);
        }
        printf("\n");
        for (int j = 0; j < number_of_Sequences; ++j) {
            if (most_similar_group_two == Aligning_Sequence[j].label) printf("Sequence-%d  ", j);
        }
        printf("\n");
        AlignAndUpdate(most_similar_group_one, most_similar_group_two);
    }
    return 0;
}

vector<string> MSA::output()
{
    vector<string> temp;
    for (int i = 0; i < number_of_Sequences; ++i) {
        temp.push_back(Aligning_Sequence[i].sequence);
    }
    return temp;
}