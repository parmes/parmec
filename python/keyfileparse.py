""" Parser for LS-DYNA GCORE keyfiles.

    NB: This only supports fixed-format keyfiles. Free- or long-format keyfiles cannot be read.
    
    To use from the command line run:
            python keyfileparse.py [-p] path_to_keyfile
    
    where:
        - if the optional "-p" switch is present then timing profile information will be output.
        - path_to_keyfile is the path to the keyfile.
    
    Once the keyfile has been parsed (this may take a while) you will have an interactive python session,
    with the variable 'keyfile' being a Keyfile object (see documentation for this below).
    
    Equivalent usage in a script is:
    
        from keyfileparse import Keyfile
        keyfile = Keyfile(path_to_keyfile)
    
    ***************************************************************************
    *  See the docstrings for the Keyfile and Card objects for the interface. *
    ***************************************************************************
    
    If the parser does not understand a particular keyfile card it will be
    skipped and a message printed.
    
    To extend the keyfile card formats the parser understands you should just be able
    to add to CARDDEFINITIONS below.
    
"""

__version__ = '1.0'

import itertools, sys, code # stdlib
import cProfile

# --- Keyfile card formats --------------------------------------------------------------------------------

# The LS-DYNA keywords which can be read by the parser is determined by the format definitions given
# in CARDDEFINITIONS below - this won't need modifying unless you need to read additional keywords.
# 
# Each keyword:format pair specifies a single data block for that keyword.
# Each line in the format specifies the fields for a single card in the data block, space-separated.
# Each field is defined in the format "N-TL" where:
#   N is the field name
#   T is the field type (see constant CONVERTERS)
#   L (optional) is the field length, otherwise set to constant DEFAULT_FIELD_LEN
#
# A *field* can be marked as optional by surrounding it's definition by brackets, e.g. "(N-TL)". If
# no fields are optional the number of fields defined must match the number found.
# If a *single* data block is itself variable length, e.g. as for DEFINE_CURVE, this is indicated by
# the final card being "..." . Note this is *not* the same as repeating datablocks under one
# keyword - e.g. as the Generator happens to do for *NODE - this is possible for any keyword.
# If a data block is variable-length then it cannot be repeated - the block is ended by the next keyword
# (this is as per the LS-DYNA manual).
#
# ***IMPORTANT ***
# The format definitions given here have been modified to match what the 4x4x4.key actually has,
# rather than strictly adhering to the manual definitions. In particular:
#   - all A8's have been changed to use the default length of 10

CARDDEFINITIONS = {'TITLE':"""TITLE-A80""",
        
                 'PART':"""HEADING-A80
                           PID-I SECID-A MID-A""",
                 
                 'PART_INERTIA':"""HEADING-A80
                                   PID-I SECID-A MID-A
                                   XC-F  YC-F  ZC-F  TM-F  IRCS-I NODEID-I
                                   IXX-F IXY-F IXZ-F IYY-F IYZ-F  IZZ-F
                                   VTX-F VTY-F VTZ-F VRX-F VRY-F VRZ-F
                                   (XL-F) (YL-F) (ZL-F) (XLIP-F) (YLIP-F) (ZLIP-F) CID-I""",
                                   # note the 1st 6 fields on the last card are only strictly option if CID is used
        
                 'CONSTRAINED_EXTRA_NODES_SET':"""PID-I NSID-I""",
        
                 'MAT_RIGID':"""MID-A10 RO-F E-F PR-F N-F COUPLE-F M-F
                                CMO-F CON1-F CON2-F
                                LCO-F A2-F A3-F V1-F V2-F V3-F""",
        
                 'SECTION_BEAM':"""SECID-A1 ELFORM-I SHRF-F QR/IRID-F CST-F SCOOR-F
                                   A-F ISS-F ITT-F J-F SA-F""",
        
                 'SECTION_DISCRETE':"""SECID-A1 DRO-I KD-F V0-F CL-F FD-F
                                       CDL-F TDL-F""",

                 'DEFINE_COORDINATE_SYSTEM':"""CID-I X0-F Y0-F Z0-F XL-F YL-F ZL-F (CIDL-I)
                                               XP-F YP-F ZP-F""",

                 'BOUNDARY_PRESCRIBED_MOTION_NODE':"""NID-I DOF-I VAD-I LCID-I""",

                 'DEFINE_SD_ORIENTATION':"""VID-I IOP-I XT-F YT-F ZT-F""",
                                       
                 'SET_NODE_LIST':"""SID-I
                                    NID1-I (NID2-I) (NID3-I) (NID4-I) (NID5-I) (NID6-I) (NID7-I) (NID8-I)
                                    ...""",
                                    
                 'NODE':"""NID-I8 X-F16 Y-F16 Z-F16 (TC-I8) (RC-I8)""", # in manual are all -10,  TC and RC are F
        
                 'DEFINE_CURVE':"""LCID-I
                                   A1-F20 O1-F20
                                   ...""",
        
                 'MAT_SPRING_NONLINEAR_ELASTIC':"""MID-A LCD-I LCR-I""", #MID A8 in manual

                 'DAMPING_GLOBAL':"""LCID-I VALDMP-F STX-F STY-F STZ-F SRX-F SRY-F SRZ-F""",

                 'LOAD_BODY_X':"""LCID-I""",

                 'LOAD_BODY_Y':"""LCID-I""",

                 'LOAD_BODY_Z':"""LCID-I""",

                 'CONTROL_TERMINATION':"""ENDTIM-F ENDCYC-I DTMIN-F ENDENG-F ENDMAS-F""",

                 'CONTROL_TIMESTEP':"""DTINIT-F TSSFAC-F ISDO-I TSLIM-F DT2MS-F LCTM-I ERODE-I MS1ST-I""",
                 
                 'ELEMENT_DISCRETE':"""EID-I8 PID-I8 N1-I8 N2-I8 VID-I8 S-F16 PF-I8 OFFSET-F16""", # All default field len in manual n10 rather than n8)
        }
        
DEFAULT_FIELD_LEN = 10

CONVERTERS = {'I':int, 'F':float, 'A':lambda a: a.strip()}

class Keyfile(object):
    """ Access an LS-DYNA **fixed-format** keyfile.
        
        This class provides access to Card objects defining the information parsed from
        the keyfile.
        
        To parse the keyfile create a Keyfile object by providing the path:
        
            keyfile = Keyfile(path)
        
        The Cards can then be accessed in a number of ways:
            
            1. Return a list of all cards, in the order they were read from the keyfile:
            
                card_list = keyfile.all_cards
            
            2. Return a dict of all cards:
            
                card_dict = keyfile.cards
                
                The keys to the dict are unique card types (i.e. keywords) and each value is a list of Card
                objects having that keyword.
                
            3. Return a single card which has the specified keyword and also all specified fieldname/value pairs
               (or None if no such card is found):
            
                card = key.getcard(keyword, field1=value1, field2=value2, ...)
            
            2. Return a list of cards with the specified (case-insenstive) keyword:
                
                cards = keyfile['keyword']
        
        Some commonly-used attributes can also be accessed using the following dictionaries:
            
            keyfile.NODES[nid]
                key is nid (int), value is tuple of floats (cx, cy, cz)
                
            self.CONSTRAINED_EXTRA_NODES_SETS
                key is pid (int), value is nsid (int)

            self.SET_NODE_LISTS
                key is sid (int), value is flat tuple of all nids
                
        Note these attributes are cached and hence much faster than using the corresponding keyfile.getcard() calls.
        
    """
    
    def __init__(self, path):

        self.cards = {} # key=keyword, value=list of cards
        self.all_cards = [] # list of all cards
        
        currblock = ''
        currkeyword = ''
        currblocklineno = None
        
        self.cardformats = self.parseformats(CARDDEFINITIONS)
        
        with open(path, 'r') as keyfile:
            for ln, line in enumerate(keyfile):
                
                if line.startswith('$'):
                    pass # skip comment lines - NB blank lines are important, they are blank cards!
                
                elif line.startswith('*') or line == '': # '' is EOF
                    
                    if currblock: # i.e. only if current block is not empty...
                        # parse it:
                        newcards = self.newcard(currkeyword, currblock, currblocklineno)
                        
                        # add the cards by keyword:
                        if currkeyword in self.cards:
                            self.cards[currkeyword].extend(newcards)
                        else:
                            self.cards[currkeyword] = newcards
                        
                        # add it in order:
                        self.all_cards.extend(newcards)
                        
                    # start a new block:
                    currblocklineno = ln + 1 # +1 as iter starts at 0
                    currblock = ''
                    currkeyword = line.lstrip('*').strip()
                
                else:
                    currblock += line
            
            # create cached items
            self.cache()
    
    def parseformats(self, carddefstr):
        """ Return a datastructure with card formats:
                key=keyword, value=(formatlines, is_variable_length, nformatlines)
        """
        
        cardformats = {}
        
        for keyword in carddefstr:
            
            # remove leading/trailing newlines which may be present due to """-quoted strings
            formatlines = carddefstr[keyword].strip('\n')
            # break the format description into lines:
            formatlines = formatlines.split('\n')
        
            if '...' in formatlines[-1]:
                is_variable_length = True
                del formatlines[-1]
            else:
                is_variable_length = False
            nformatlines = len(formatlines)
            
            newformat = (formatlines, is_variable_length, nformatlines)

            cardformats[keyword] = newformat

        return cardformats
    
    def getcard(self, keyword, **kwargs):
        """ Return a single card with the specified keyword and matching the fieldname=value pairs
            given as kwargs.
            Keys are case-insensitive
        """
        
        keyword = keyword.upper()
        if keyword not in self.cards:
            return None
        
        kwargs_set = set((k.upper(), v) for (k, v) in kwargs.items())
        
        for c in self.cards[keyword]:
            c_set = set(c.fields.items())
            if kwargs_set.issubset(c_set):
                return c
        return None

    def __getitem__(self, keyword):
        """ return a list of Cards having this keyword (case-insensitive) """
        try:
            return self.cards[keyword.upper()]
        except KeyError:
            raise KeyError('no cards for keyword "%s" found' % (self.keyword))

    def cache(self):
        """ Provide cached attributes as defined below for speed:
        
            self.NODES                             # key=nid, value=(cx, cy, cz)
            self.CONSTRAINED_EXTRA_NODES_SETS = {} # key=pid, value=nsid    
            self.SET_NODE_LISTS = {}               # key=sid, value = flat tuple of all nids
        """
        
        self.NODES = {}
        for card in self['NODE']:
            self.NODES[card['nid']] = (card['x'], card['y'], card['z'])
            
        self.CONSTRAINED_EXTRA_NODES_SETS = {}
        for card in self['CONSTRAINED_EXTRA_NODES_SET']:
            self.CONSTRAINED_EXTRA_NODES_SETS[card['pid']] = card['nsid']
        
        self.SET_NODE_LISTS = {}
        for card in self['SET_NODE_LIST']:
            nids = flatten([v for (f, v) in card.fields.items() if f.startswith('NID')])
            self.SET_NODE_LISTS[card['sid']] = nids

    def __str__(self):
        output = []
        for kw in sorted(self.cards):
            output.extend([str(c) for c in self.cards[kw]])
        return '\n\n'.join(output)

    def newcard(self, keyword, datalines, lineno=None):
        """ Create new Card objects from keyfile data.
        
            keyword: str
            datalines: multi-line str
            lineno: the line number in the file for keyword (not required but improves error outputs)
            
            Returns a list of Card objects, potentially empty if none were generated.
        """
        
        cards = []
        
        # check keyword is defined:
        if keyword not in self.cardformats:
            print 'skipping %s' % keyword
            return cards
            
        (formatlines, is_variable_length, nformatlines) = self.cardformats[keyword]
        
        # prepare data lines:
        datalines = datalines.rstrip('\n') # only remove the final newline on the block
        # CAUTION: this is not totally safe as \n\n would also be removed, and the final card could be blank
        # but I have never seen this
        datalines = datalines.split('\n')
        nlines = len(datalines)
        if nlines == nformatlines or is_variable_length: # keyword has a single block of data lines
            blocks = [datalines]
        else: # keyword has repeated blocks
            if not (nlines % nformatlines) == 0: # check have whole repeats
                print 'parsing blocks for keyword "%s" at line %i failed:' % (keyword, lineno)
                print '%i lines in block but %i lines in format definition - not compatible' % (nlines, nformatlines)
                raise ValueError()
                
            # group lines them into format-length blocks:
            blocks = zip(*[iter(datalines)]*nformatlines) # see docs for zip() for this!
            assert len(blocks) == nlines / nformatlines
        
        # now process each block
        for block in blocks:
            #print 'new block:', block
            fields = {}
            
            # iterate over data lines and format description lines in parallel:
            for ln, dataline in enumerate(block):
                
                if is_variable_length and ln > (nformatlines-1):
                    linefmt = formatlines[-1]   # use last format line repeatedly
                else:
                    linefmt = formatlines[ln]
                                        
                # iterate across fields and fieldcodes:
                for fieldcode in linefmt.split():
                    
                    # is field optional, e.g. "(NID2-I)":
                    if fieldcode.startswith('(') and fieldcode.endswith(')'):
                        is_optional = True
                        fieldcode = fieldcode.strip('()')
                    else:
                        is_optional = False
                    
                    fname, ftype = fieldcode.split('-')
                    fname = fname.upper()
                    
                    # get field length:
                    if len(ftype) > 1:
                        flen = int(ftype[1:])
                    else:
                        flen = DEFAULT_FIELD_LEN
                    ftype = ftype[0]
                    fieldstr = dataline[0:flen] # NB if flen is off the end of the line then this will just go to the end
                    dataline = dataline[flen:]
                    
                    # check whether there is data to convert:
                    if is_optional and (fieldstr.isspace() or fieldstr == ''):
                        continue # this field isn't specified
                    
                    # convert the field:
                    converter = CONVERTERS[ftype]
                
                    try:
                        val = converter(fieldstr)
                    except Exception as e:
                        print 'conversion of field "%s" in keyword "%s" at line %i failed:' % (fname, keyword, lineno)
                        print 'fields read ok:%s' % fields
                        print 'fieldstring:%r, fieldname:%r, type:%r, len:%r, is_optional:%s' % (fieldstr, fname, ftype, flen, is_optional)
                        print 'remaining line: %r' % dataline
                        raise e
                    
                    if is_variable_length and ln == (nformatlines-1): # on last format line:
                        fields[fname] = [val]
                    elif is_variable_length and ln > (nformatlines-1): # on repeated lines:
                        fields[fname].append(val)
                    else:   # non-variable line in variable block, or non-variable block:
                        fields[fname] = val
    
                
            # now convert any list-fields to tuples (required so that values are hashable)
            for f, v in fields.items():
                if isinstance(v, list):
                    fields[f] = tuple(v)
                    
            cards.append(Card(keyword, fields))
            
        return cards

def flatten(seq_of_seqs):
    return list(itertools.chain(*seq_of_seqs))
    
class Card(object):
    """ Represents a single keyfile card.
    
        Access to individial fields is via dictionary-style access using the (case-insensitive) fieldname, e.g. 
        
            card['NODE'] or card['node']
        
        The dictionary of all fields is available as
            card.fields
        
        The keyword for the card is provided as field "KEYWORD".
    """
    
    def __init__(self, keyword, fields_dict):
        self.keyword = keyword.upper()
        self.fields = fields_dict
        self.fields['KEYWORD'] = self.keyword
    
    def __getitem__(self, fieldname):
        """ Allow access to field as card[fieldname] (case-insensitive) """
        try:
            return self.fields[fieldname.upper()]
        except:
            raise KeyError('%s has no field "%s"' % (self.keyword, fieldname.upper()))
            
    def __str__(self):
        return '*%s:%s' % (self.keyword, self.fields)
    
if __name__ == "__main__":
    if sys.argv[1] == '-p':
        keypath = sys.argv[2]
        cProfile.run("keyfile = Keyfile(keypath)", sort='cumulative')
    else:
        keypath = sys.argv[1]
        keyfile = Keyfile(keypath)
    code.interact('keyfile loaded as variable "keyfile"', local=locals())
    
