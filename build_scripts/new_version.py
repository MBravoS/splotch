"""


Author: Robin Cook, Matias Bravo
Date: 09/08/19
"""

#from sys import argv

#version_update=int(argv[1])
#changelog_message=argv[2]

changelog = "../changelog.txt"
readme = "../README.md"
setup = "../setup.py"
conf = "../docs/conf.py"

import os
from argparse import ArgumentParser, RawTextHelpFormatter


def get_current():
	# Returns current version
	pass
	


parser = ArgumentParser(description="Parse argument for changing versions:\n" 
									"  * Version (+v/--version): Reserved for largest updates when moving to new versions.\n"
									"  * Major (+M/--major): When a major update/new feature has been made to the code.\n"
									"  * Minor (+m/--minor): The addition of small changes to an existing function.\n"
									"  * Bug (+b/--bug): Small adjustments and bug fixes.\n", formatter_class=RawTextHelpFormatter, prefix_chars='-+')
parser.add_argument('-f','--full',
					action='store',dest='full',type=str,default=False,
					help="Change the full version number X.X.X.X",metavar="VERSION")
parser.add_argument('+v','--version',
					action='store',dest='version',nargs='?',type=int,default=False,const=True,
					help="Increase the version number X.0.0.0")
parser.add_argument('+M','--major',
					action='store',dest='major',nargs='?',type=int,default=False,const=True,
					help="Increase the major version number -.X.0.0",metavar="MAJOR")
parser.add_argument('+m','--minor',
					action='store',dest='minor',nargs='?',type=int,default=False,const=True,
					help="Increase the minor version number -.-.X.0",metavar="MINOR")
parser.add_argument('+b','--bug',
					action='store',dest='bug',nargs='?',type=int,default=False,const=True,
					help="Increase the minor version number -.-.-.X",metavar="BUG")
parser.add_argument('-m','--message',
					action='store',dest='message',type=str,default=None,
					help="The message to be added to changelog.txt. If no --message given, do not add to changelog.txt",metavar="MESSAGE")
parser.add_argument('-c','--current',
					action='store_true',dest='current',default=False,
					help="Print the current version and exit.")

args = parser.parse_args()

# Get the current version
currVersion = [2,4,5,7]#get_current()
currVersionStr = ".".join([str(kk) for kk in currVersion])
print(f"Current version is: {currVersionStr}")

# Print current version and YEET on out of there
if (args.current == True):
	print("The current version is {0}".format('ding'))
	exit()

# Validate inputs
nArgs = sum([not args.full is False, not args.version is False, not args.major is False, not args.minor is False, not args.bug is False ])

print(args.version)

if (nArgs == 1):
	if (args.full):
		newVersion = args.full.split('.')
		if (len(newVersion) != 4 or len(args.full) < 7):
			print(f"\nERROR: --full argument not valid, must have the format X.X.X.X, but {args.full} was given.\n  -- ABORTING --\n"); exit()
	elif (not args.version is False):
		newVersion = [currVersion[0]+args.version,0,0,0] if (args.version is True) else [args.version,0,0,0]
	elif (not args.major is False):
		newVersion = [currVersion[0],currVersion[1]+args.major,0,0] if (args.major is True) else [currVersion[0],args.major,0,0]
	elif (not args.minor is False):
		newVersion = [currVersion[0],currVersion[1],currVersion[2]+args.minor,0] if (args.minor is True) else [currVersion[0],currVersion[1],args.minor,0]
	elif (not args.bug is False): 
		newVersion = [currVersion[0],currVersion[1],currVersion[2],currVersion[3]+args.bug] if (args.bug is True) else [currVersion[0],currVersion[1],currVersion[2],args.bug]

elif (nArgs == 0):
	print("\nERROR: Nothing specified to change version number. Use one of [+v.+M.+m.+b] | -f.\n  -- ABORTING --\n"); exit()
else:
	print("\nERROR: Too many arguments given. Use one of [+v.+M.+m.+b] | -f.\n  -- ABORTING --\n"); exit()


newVersionStr = ".".join([str(kk) for kk in newVersion])
print(f"\nNew version is {newVersionStr}")