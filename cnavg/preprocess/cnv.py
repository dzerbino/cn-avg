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

"""CNV data holder"""

import sys
import cnavg.basics.coords as coords
from breakend import Breakend

################################################
## CNV
################################################

class CNV(coords.Region):
	""" CNV data holder"""
        def __init__(self, chr, start, finish, val, name, softStart=None, softFinish=None):
		super(CNV, self).__init__(chr, start, finish)
                self.val = list(val)
		self.softStart = softStart
		self.softFinish = softFinish
		self.startCap = None
		self.finishCap = None
		self.name = str(name)

        def __str__(self):
		if self.softStart is None or self.softFinish is None:
			return "[%s] %s:%i-%i=%s" % (self.name, self.chr, self.start, self.finish, str(self.val))
		else:
			return "[%s] %s:%i/%li-%i/%li=%s" % (self.name, self.chr, self.start, self.softStart, self.softFinish, self.finish, str(self.val))

	def addStartCap(self, graph):
		self.startCap = Breakend(self.chr, self.start, True, self.name + ":start", finish=self.softStart)
		self.startCap.createPartner(graph)
		assert self.val >= 0
		self.startCap.segment = self.val
		graph.addBreakend(self.startCap)

	def addFinishCap(self, graph):
		if self.softFinish is not None:
			self.finishCap = Breakend(self.chr, self.softFinish, False, self.name + ":finish", finish=self.finish)
		else:
			self.finishCap = Breakend(self.chr, self.finish - 1, False, self.name + ":finish", finish=self.finish)
		self.finishCap.createPartner(graph)
		assert self.val >= 0
		self.finishCap.segment = self.val
		graph.addBreakend(self.finishCap)

	def validate(self):
		super(CNV, self).validate()
		if self.softStart is not None:
			assert self.softStart >= self.start
			assert self.softStart <= self.finish
		if self.softFinish is not None:
			assert self.softFinish >= self.start
			assert self.softFinish <= self.finish
		if self.softStart is not None and self.softFinish is not None:
			assert self.softFinish >= self.softStart
		return True

	def ploidy(self):
		return len(self.val)
