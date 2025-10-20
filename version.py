# convert git shortehand hash into version string
import subprocess

git = ["git",  "log",  "--pretty=format:%h", "-n", "1"]
p = subprocess.Popen(git, stdout=subprocess.PIPE)
commit = p.communicate()[0].decode('utf-8').strip()
git = ["git", "show", "-s", "--format=%cd", "--date=short", commit]
p = subprocess.Popen(git, stdout=subprocess.PIPE)
date = p.communicate()[0].decode('utf-8').replace('\n','')

# Read existing version.h if it exists
comm0 = ""
try:
    with open('version.h') as f:
        lines = f.readlines()
        if lines:
            # Extract existing hash, handling both quoted and unquoted formats
            first_line = lines[0].strip()
            if '#define VERSION_HASH' in first_line:
                comm0 = first_line.split()[-1].strip().strip("'\"")
except FileNotFoundError:
    pass  # It's OK if the file doesn't exist yet

# Update the version file if commit has changed or if format needs fixing
if comm0 != commit:
    with open('version.h','w') as f:
        f.write('#define VERSION_HASH "%s"\n' % commit)
        f.write('#define VERSION_DATE "%s"\n' % date)