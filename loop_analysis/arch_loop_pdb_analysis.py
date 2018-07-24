import pyrosetta
from pyrosetta.rosetta import *
from Blueprint import Blueprint
import numpy as np
import sys,re,glob,os

pyrosetta.init()

def MakeBlueprint(pose):
	DSSP=protocols.moves.DsspMover()
        DSSP.apply(pose)  
        ss = pose.secstruct()
	seq = pose.sequence()
        # Get abegos
        abego = core.sequence.get_abego(pose) 
        abego = list(abego)     
        res_list=[]
        for k,ss_res in enumerate(ss):
                res_list.append( [k+1,seq[k],ss_res+abego[k],'.'] )

        blue = Blueprint(data=res_list)
        return blue

def Angle(v1,v2):
    scalar = v1.dot(v2) / ( np.linalg.norm(v1) * np.linalg.norm(v2) )
    return np.rad2deg( np.arccos(scalar) )

def StrandDirectionTriad(pose,triad):
        a,b,c = triad
        cen1 = np.array( pose.residue( a ).xyz("CA") + pose.residue( b ).xyz("CA") ) / 2.0
        cen2 = np.array( pose.residue( b ).xyz("CA") + pose.residue( c ).xyz("CA") ) / 2.0
        v = cen2-cen1
        return v


def StrandPairing(pose):
    blue = MakeBlueprint(pose)
    dssp=core.scoring.dssp.Dssp(pose)
    pairset=core.scoring.dssp.StrandPairingSet(pose)
    strand_dic_pairs={} ; strand_dic_a_pairs={} ; strand_dic_p_pairs={} ; sspair=[]
    for i in range( 1, pairset.size()+1 ):
        sp=pairset.strand_pairing(i)
        start=sp.begin1()
        end=sp.end1()
        segment = blue.residue_segment(start)
        partner = sp.get_pair(start)
        partner2 = sp.get_pair(end)
        segment2 = blue.residue_segment(partner)
        sspair.append([segment,segment2])
    return sspair

def AbegoConnectionType(pdbname, connection):
	pose=core.import_pose.pose_from_file(pdbname)

	polyala = protocols.pose_creation.MakePolyXMover(aa="ALA",keep_pro=0,keep_gly=0,keep_disulfide_cys=0)
	polyala.apply(pose)

	nres = pose.total_residue()

	# Make blueprint
	blue = MakeBlueprint(pose)

	r = re.compile('([HEL]\d+)-?')
	seg_list = r.findall(blue.topology())
        r = re.compile('([E]\d+)-?')
        strand_list = r.findall(blue.topology())	

        # Strand pairing list
	sspair_lst = StrandPairing(pose)

	st=''
	for strand1,strand2 in zip(strand_list[:-1],strand_list[1:]):
		seg1 = blue.segment_dict[strand1]
		seg2 = blue.segment_dict[strand2]
        	strand1_idx = seg_list.index(strand1)
	        strand2_idx = seg_list.index(strand2) 
		loop_data=[]
                if ( len(seg1.bp_data)<4 ) or ( len(seg2.bp_data)<4 ):
                        continue

		if ( strand2_idx - strand1_idx  == 2 ):
	                for idx in range(strand1_idx+1,strand2_idx):
        	            segment = seg_list[idx]
                	    seg = blue.segment_dict[segment]
	                    loop_data.extend( seg.bp_data )

		elif ( strand2_idx - strand1_idx  > 2 ):
			if abs( ( seg1.bp_data[-1][0]+1 ) - (seg2.bp_data[0][0]-1 ) ) + 1 <=6:
				loop_data =[]
		                for idx in range(strand1_idx+1,strand2_idx):
                		    segment = seg_list[idx]
		                    seg = blue.segment_dict[segment]
                		    loop_data.extend( seg.bp_data )

		if len(loop_data) > 0:				
				abego='' ; seq=''
				first_residue_index = loop_data[0][0] - 1
				seq2 = blue.bp_data[first_residue_index-1][1][0] # ncap
				abego2=blue.bp_data[first_residue_index-1][2][-1]
				for res in loop_data:
					abego += res[2][-1]
					seq += res[1][0]
					pos = res[0]-1
					seq2 += res[1][0]
					abego2+=res[2][-1]

				seq2 += blue.bp_data[pos+1][1][0]
				abego2 += blue.bp_data[pos+1][2][-1]

				# Check whether EE loop is hairpin or arch
				if connection == 'EE':
					loop_type='arch'
			                for sspair in sspair_lst:
                        			if strand1 in sspair and strand2 in sspair:
			                            loop_type='hairpin'

					# compute angle of strand directionality
					triad = [loop_data[0][0]-3,loop_data[0][0]-2,loop_data[0][0]-1]
                                        v1 = StrandDirectionTriad(pose,triad)
                                        triad = [loop_data[-1][0]+1,loop_data[-1][0]+2,loop_data[-1][0]+3]
                                        v2 = StrandDirectionTriad(pose,triad)

					trans = pose.residue(loop_data[-1][0]+1).xyz("CA")- pose.residue(loop_data[0][0]-1).xyz("CA")
					# check concavity.
					conc_v = -v1+v2
					cacb1 = pose.residue(loop_data[0][0]-1).xyz("CB")- pose.residue(loop_data[0][0]-1).xyz("CA")
                                        cacb2 = pose.residue(loop_data[-1][0]+1).xyz("CB")- pose.residue(loop_data[-1][0]+1).xyz("CA")
					conc1=False ; conc2=False
					if cacb1.dot(trans) > 0:
						conc1 = True
                                        if cacb2.dot(trans.negate()) > 0:
                                                conc2 = True	

					# loop center of mass
					cm=pose.residue(loop_data[0][0]).xyz("CA")  
					if len(loop_data)>1:
						for res in loop_data[1:]:
							cm =  cm + pose.residue(res[0]).xyz("CA")
					cm /= len(loop_data)

					# previous position
					pos1 = np.array( pose.residue(loop_data[0][0]-1).xyz("CA") + pose.residue(loop_data[0][0]-2).xyz("CA") ) / 2.0 
					# next position
					pos2 = np.array( pose.residue(loop_data[-1][0]+1).xyz("CA") + pose.residue(loop_data[-1][0]+2).xyz("CA") ) / 2.0
					va = pos1-cm
					vb = pos2-cm
					bend = Angle(va,vb)
							

					# Calculate geometry of connection
					angle = Angle(v1,v2)
					dist = pose.residue(loop_data[0][0]-1).xyz("CA").distance( pose.residue(loop_data[-1][0]+1).xyz("CA")  )

					# calculate dihedral
					pos1 = pose.residue(loop_data[0][0]-3).xyz("CA")
					pos2 = pose.residue(loop_data[0][0]-1).xyz("CA")
					pos3 = pose.residue(loop_data[-1][0]+1).xyz("CA")
					pos4 = pose.residue(loop_data[-1][0]+3).xyz("CA")
					dih=numeric.dihedral_degrees(pos1,pos2,pos3,pos4)

					st+='%s %s %s %s %s %s %s %s %s %.3f %.3f %.3f %.3f %s %s\n' %(pdbname, connection, loop_data[0][0], loop_data[-1][0], seq, seq2, abego, abego2, loop_type, dist, angle, bend, dih, conc1, conc2)

	return st

##############################
filein = open( sys.argv[1] ) # list of pdbfiles to analyze
fileout = open( sys.argv[1]+'-OUT.log','w')
for line in filein:
	pdbname = line.split()[0]
	st=AbegoConnectionType(pdbname, 'EE')
	fileout.write(st)

fileout.close()
##############################
