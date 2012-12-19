# Copyright (c) 2012, Daniel Zerbino
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# (1) Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer. 
# 
# (2) Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in
# the documentation and/or other materials provided with the
# distribution.  
# 
# (3)The name of the author may not be used to
# endorse or promote products derived from this software without
# specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#!/usr/bin/env python

"""
Definition of geometric objects along chromosomes.
"""

class Region(object):
       """ 
       A contiguous region.
       """
       def __init__(self, chr, start, finish):
              self.chr = str(chr)
              self.start = int(start)
              self.finish = int(finish)

       def __str__(self):
              return "\t".join(map(str, [self.chr, self.start, self.finish]))

       def __cmp__(self, other):
              if cmp(self.chr, other.chr):
                     return cmp(self.chr, other.chr)
              elif self.start > other.finish:
                     return 1
              elif self.finish < other.start:
                     return -1
              else:
                     return 0

       def validate(self):
              """ Validation """
              assert self.finish >= self.start
              return True

       def length(self):
              return self.finish - self.start + 1

###############################################
## Position
###############################################

class Position(Region):
       """ A 1-bp region """
       def __init__(self, chr, pos):
              super(Position, self).__init__(chr, pos, pos)
              self.pos = int(pos)

###############################################
## Oriented Region
###############################################

def orientationChar(orientation):
       if orientation:
              return '+'
       else:
              return '-'

class OrientedRegion(Region):
       """ A region with an assigned strand """
       def __init__(self, chr, start, finish, orientation):
              super(OrientedRegion, self).__init__(chr, start, finish)
              self.orientation = bool(orientation)

       def position(self):
              """ Left-most coordinate """
              if self.orientation:
                     return self.start
              else:
                     return self.finish

       def __cmp__(self, other):
              res = super(OrientedRegion, self).__cmp__(other) 
              if res != 0:
                     return res
              else:
                     return cmp(self.orientation, other.orientation)

       def __str__(self):
              return "\t".join(map(str, [self.chr, self.start, self.finish, orientationChar(self.orientation)]))

###############################################
## Oriented Position
###############################################

class OrientedPosition(OrientedRegion):
       """ A position with an assigned strand """
       def __init__(self, chr, pos, orientation):
              super(OrientedPosition, self).__init__(chr, pos, pos, orientation)
              self.pos = int(pos)
