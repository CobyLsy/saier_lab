#!/usr/bin/env python

'''
#Test Cases:
#q: 1.A.95.2.5-YP_009361958
#t: 1.A.100.1.1-B2X7D9

#given by mmseqs
qaln = 'MGFDIGGDIGKPLKDAFDKFGADIKMTFLTVLNWMK--WISIG------ILIVISVI-------LICKIIKVLFQCGKCLLSCFGFCKK'
taln = 'MGFSINFD---PIINKFREFQTNINHNINEQLDKLKMVWINLGSHIKYWFIIIISILTILFILFLLIKITKLILNCKKIFSCCCNVCCK'

#given by mmseqs
qlen = 120
tlen = 127

#rawqTMS = [[26, 49], [54,72]] #1
#rawqTMS = [[6,15], [36,42], [47,57]] #3
#rawqTMS = [[13,31], [39,57]] #4
#rawqTMS =[[1,20],[71,80]] #5
#rawqTMS = [[12,35], [46,54]] #6
rawqTMS = [[8,30], [36,50], [59,73]]

#rawtTMS = [[68, 92]] #1
#rawtTMS = [[27,36], [52,90]] #3
#rawtTMS = [[26,38], [42,75], [83,104]] #4
#rawtTMS =[[14,29],[102,111]] #5
#rawtTMS = [[26,39], [45,75]] #6
rawtTMS = [[24,35], [42,68], [73,96]]

#Given by mmseqs
#Already 0-indexed, but are GIVEN BY MMSEQS AS 1-INDEXED
qalnStartPos = 0
qalnEndPos = 73

talnStartPos = 21 
talnEndPos = 106 

'''

#query : WP_002485759.1
#target : 3.A.3.6.1-P20021
qaln = 'QAPQEVKKDKNIYRVEGFSCANCAGKFERNVKKIPGVEDAKVNFGASKISVYGEATIEELEKAGAFENLKVASEKPVRQATQEINQEKEDEKEEKVPFYKKHSTLLYSTLLIVFGYLSVFVNGDENIVTTLLFVASMLIGGLSLFKVGFQNLLRFEFDMKTLMTVAVIGGAIIGEWAEVSVVVILFAISEALERFSMDKARQSIRSLMDIAPKEALVRRKGQEMMVHVDDIAVGDIMIVKPGQKIAMDGMVVSGYSAVNQAAITGESVPVEKAVDDEVFAGTLNEEGLLEVEITKLVEDTTISKIIHLVEEAQGERAPSQAFVEKFAKYYTPIIMIIAALVAVVPPLFFGASWETWVYQGLAVLVVGCPCALVISTPISIVSAIGNAAKKGVLIKGGVYLEEMGALKAIAFDKTGTLTKGVPVVTDFNVLNKQVDENEMLSIITALEYRSQHPLASAIMKRAEEANISYSDVVIDDFSSITGKGIKGTVDGTTYYIGSPKLFKELSNSSFDKNLEKKVATLQNQGKTAMVVGTDKEILAIIAVADEVRESSKEVIQKLHQLGIKNTIMLTGDNKGTANAIGSHVGVKEVQAELMPQDKLDYIKQLKSEYNNVAMIGDGVNDAPALAASTVGIAMGGAGTDTALETADVALMGDDLRKLPFTVKLSRKALNIIKANITFAIAIKFIALLLVIPGWLTLWIAILSDMGATLLVALNSLRLMRVKD'
taln = 'QKVKLMEEEMNVYRVQGFTCANCAGKFEKNVKKIPGVQDAKVNFGASKIDVYGNASVEELEKAGAFENLKVSPEKLANQTIQRVKDDTKAHKEEKTPFYKKHSTLLFATLLIAFGYLSHFVNGEDNLVTSMLFVGSIVIGGYSLFKVGFQNLIRFDFDMKTLMTVAVIGATIIGKWAEASIVVILFAISEALERFSMDRSRQSIRSLMDIAPKEALVRRNGQEIIIHVDDIAVGDIMIVKPGEKIAMDGIIVNGLSAVNQAAITGESVPVSKAVDDEVFAGTLNEEGLIEVKITKYVEDTTITKIIHLVEEAQGERAPAQAFVDKFAKYYTPIIMVIAALVAVVPPLFFGGSWDTWVYQGLAVLVVGCPCALVISTPISIVSAIGNAAKKGVLVKGGVYLEKLGAIKTVAFDKTGTLTKGVPVVTDFEVLNDQVEEKELFSIITALEYRSQHPLASAIMKKAEQDNIPYSNVQVEEFTSITGRGIKGIVNGTTYYIGSPKLFKELNVSDFSLGFENNVKILQNQGKTAMIIGTEKTILGVIAVADEVRETSKNVIQKLHQLGIKQTIMLTGDNQGTANAIGTHVGVSDIQSELMPQDKLDYIKKMQSEYDNVAMIGDGVNDAPALAASTVGIAMGGAGTDTAIETADIALMGDDLSKLPFAVRLSRKTLNIIKANITFAIGIKIIALLLVIPGWLTLWIAILSDMGATILVALNSLRLMRVKD'

rawqTMS = [[184,201], [206, 223], [236, 253], [258, 275], [413, 430], [439, 463], [753, 770], [775, 792]]
rawtTMS = [[107, 124], [131, 148], [165, 189], [334, 351], [364, 388], [395, 414], [627, 646], [677, 694], [699, 716]]

qalnStartPos = 79
qalnEndPos = 801

talnStartPos = 3
talnEndPos = 725

def oneToZeroIndexed(coords):
    for i in coords:
        i[0] = i[0] - 1
        i[1] = i[1] - 1
    return coords

qTMS = oneToZeroIndexed(rawqTMS)
tTMS = oneToZeroIndexed(rawtTMS)



#Minimum Residue Number for an overlap to be counted
minRes = 8

#Assume coords have been already corrected to 0-based indeces
def relRMScoordsi (coords, alnStart, alnEnd):
	#Remove TMSs that are not in the alignment region
	#using the equation:
	#   relTMCScoord = hmmtopCoord - alnStart
    returnTMS = []
    for i in coords:
        i[0] = i[0] - alnStart
        i[1] = i[1] - alnStart

        if i [1] <= 0:
            continue

        elif i[0] >= alnEnd - alnStart:
            continue
        
        else:
            returnTMS.append(i)

    return returnTMS
	#Return only the TMS that are in the alignment.

qTMSWithin = relRMScoordsi (qTMS, qalnStartPos, qalnEndPos)
tTMSWithin = relRMScoordsi (tTMS, talnStartPos, talnEndPos)

print(qTMSWithin)
print(tTMSWithin)

#Helper function that returns true if given starting index is within the range of a TMS, otherwise returns false
def findTMS(givenIndex, TMSList):
    for i in TMSList:
        if i[0] <= givenIndex and i[1] >= givenIndex:
            return True
    
    return False

#Helper function that remove a TMS from the list if given Index is within the range of a TMS
def removeTMS(givenIndex, TMSList):
    for i in TMSList:
        if i[0] <= givenIndex and i[1] >= givenIndex:
            TMSList.remove(i) #CHECK SYNTAX
            return

def overlapScore (qTMSWithin, tTMSWithin):
    overlaps = 0
    
    qTracking = 0
    tTracking = 0

    #copy two lists of TMSWithin, SYNTAX MIGHT NOT BE CORRECT!!!!
    tempQWithin = qTMSWithin
    tempTWithin = tTMSWithin

    residueCount = 0
    for i in range (0, len(qaln)): #aligned region has the same length for query and target (counting the gaps)

        #if encounter gaps, actual indices ("Tracking") don't increment, and no need to compare TMS
        print (str(i) + ': ' + qaln[i])
        print (str(i) + ': ' + taln[i])

        if qaln[i] == "-":
            qTracking = qTracking - 1
        if taln[i] == "-":
            tTracking = tTracking - 1
        

        #identical Residue Position Found
        if findTMS(qTracking, tempQWithin) == True and qaln[i] != "-" and taln[i] != "-":
            if findTMS(tTracking, tempTWithin) == True :
                print ([qTracking, tTracking])
                residueCount = residueCount + 1
            else:
                residueCount = 0 #making sure next residue count doesn't inherit previous count
        elif findTMS(qTracking, tempQWithin) == False :
            residueCount = 0 #making sure next residue count doesn't inherit previous count

        if residueCount == minRes:
            overlaps = overlaps + 1
            residueCount = 0 #making sure next residue count doesn't inherit previous count

            #If one pair TMS is already counted as overlapped, they can't be used for future comparison
            removeTMS(qTracking, tempQWithin) 
            removeTMS(tTracking, tempTWithin)

        #tracking actual indices
        qTracking = qTracking + 1
        tTracking = tTracking + 1

    #print(tempQWithin)
    #print(tempTWithin)
    return overlaps

#print(qTMSWithin)
#print(tTMSWithin)
print("Overlaps:" + str(overlapScore(qTMSWithin, tTMSWithin)))
        

        
        





