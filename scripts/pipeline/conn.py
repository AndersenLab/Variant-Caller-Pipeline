# First - check whether the files in question exist.
import os, re
from peewee import *
from datetime import datetime
import subprocess

## Fix for python 2.6

if "check_output" not in dir( subprocess ): # duck punch it in!
    def f(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = f

#==========#
# Database #
#==========#

# Determine system type:
system_type =  os.uname()[0]
md5_system = {"Linux" : "md5sum", "Darwin":"md5"} # Needed to make md5 work locally and on cluster.
# If local cd into appropriate directory
if system_type == "Darwin":
	os.chdir("../../data/fq")
	db_loc = "../"
else:
	db_loc = "/exports/people/andersenlab/dec211/"

db = SqliteDatabase(db_loc + "seq_data.db", threadlocals=True)

db.connect()

class EAV(Model):
    Entity_Group = CharField(index=True, null = True)
    Entity = CharField(index=True)
    Sub_Entity = CharField(index=True, null = True) # Used for option secondary super Entity_Groupings
    Attribute = CharField(index=True)
    Sub_Attribute = CharField(index=True, null = True) # Used for optional secondary sub Entity_Groupings
    Value = CharField(index=True)
    Tool = CharField(index=True, null=True)
    DateTime = DateTimeField()

    class Meta:
        database = db # This model uses the "seq_data.db" database.
        indexes = (
            # Specify a unique multi-column index on from/to-user.
            (('Entity', 'Attribute', 'Value',), True),
        )

db.create_tables([EAV], safe=True)


def save_eav(Entity, Attribute, Value, Entity_Group=None, Sub_Entity=None, Sub_Attribute=None, Tool=None):
	record = EAV(Entity=Entity, Attribute=Attribute, Value=Value, Entity_Group=Entity_Group, Sub_Entity=Sub_Entity, Sub_Attribute=Sub_Attribute, Tool=Tool, DateTime=datetime.now())
	try:
		record.save()
	except:
		record.update()

def save_md5(files = [], Entity = ""):
	print "in md5"
	md5 = subprocess.check_output("parallel %s ::: %s" % (md5_system[system_type], ' '.join(files)), shell=True)
	m = re.findall('MD5 \((.*)\) = (.*)', md5)
	for i in m:
		# Check File Hashes
		save_eav(Entity, "File Hash", i[1], Entity_Group = "MD5 Hash", Sub_Entity= i[0], Tool = "MD5")

