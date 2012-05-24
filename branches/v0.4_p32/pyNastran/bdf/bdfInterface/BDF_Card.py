## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import sys
import copy #???

def collapse(card):
    """doesnt work for fields==0"""
    raise Exception('broken')
    #print "cardStart = ",card
    card2 = []
    for imax,field in enumerate(reversed(card)):
        #if field != '':
        if field or str(field):
            break
        ###
    ###
    #print "imax = ",imax
    nStop = len(card)-imax
    #print "nStop = ",nStop

    for i,field in enumerate(card):
        if i==nStop:
            break
        card2.append(field)
    ###
    #print "cardEnd = ",card2
    return card2

class BDF_Card(object):
    def __init__(self,card=None,oldCardObj=None,debug=False):
        self.debug = debug
        if card:
            self.card = self.wipeEmptyFields(card)
            #self.card = card
            #print card
            #self.oldCard = oldCardObj
            self.nfields = len(self.card)
        else:
            self.oldCard = None
            self.card = None
            self.nfields = None
        ###
        #print "made card"

    def Is(self,cardName):
        """
        Returns True if the card is of type cardName
        @param self the object pointer
        @param cardName the cardName to compare against
        @retval IsACardName True/False
        """
        if card.field(0)==cardName.upper():
            return True
        return False

    def wipeEmptyFields(self,card):
        """
        Removes an trailing Nones from the card.  Also converts empty strings to None.
        @param self the object pointer
        @param card the fields on the card as a list
        @retval shortCard the card with no trailing blank fields
        """
        #print "cardA = ",card
        cardB = []
        for field in card:
            if isinstance(field,str) and field.strip()=='':
                field = None
            ###
            cardB.append(field)
        ###
        #print "cardB = ",cardB
        i = 0
        iMax = 0
        while i<len(card):
            if cardB[i] is not None:
                iMax = i
            ###
            i+=1
        ###
        #print "i=%s iMax=%s"%(i,iMax)
        #print "cardC = ",cardB[:iMax+1],'\n'
        
        return cardB[:iMax+1]

    def __repr__(self):
        """
        prints the card as a list
        @param self the object pointer
        @retval msg the string representation of the card
        """
        return str(self.card)

    def nFields(self):
        """
        gets how many fields are on the card
        @param self the object pointer
        @retval nFields the number of fields on the card
        """
        return self.nfields

    def fields(self,i=0,j=None,defaults=[],debug=False):
        """
        gets multiple fields on the card
        @param self the object pointer
        @param i the ith field on the card (following list notation)
        @param j the jth field on the card (None means till the end of the card)
        @param defaults the default value for the field (as a list) len(defaults)=i-j-1
        @param debug prints out the values at intermediate steps
        @retval the values on the ith-jth fields
        """
        if j is None:
            if self.nfields==None:
                return [None]
            j = self.nfields
            #print "jfield = ",j
        
        if defaults==[]:
            defaults=[None]*(j-i+1)
        out = []
        
        d = 0
        for n in range(i,j):
            if debug:
                print("  default = ",defaults[d])
            value = self.field(n,defaults[d])
            if debug:
                print("  n = ",n)
                print("  self.field(n) = ",self.field(n))
            out.append(value)
            d+=1
        ###
        #print "card.fields out = %s" %(out)
        return out

    def field(self,i,default=None):
        """
        gets the ith field on the card
        @param self the object pointer
        @param i the ith field on the card (following list notation)
        @param default the default value for the field
        @retval the value on the ith field
        """
        if i<self.nfields and self.card[i] is not None and self.card[i] is not '':
            #print "self.card[%s] = %s"%(i,str(self.card[i]))
            return self.card[i]
        else:
            return default
        ###
    ###

    def replaceExpression(self,fieldNew,fieldOld,replaceChar='=',replaceChar2=''):
        """used for nastran = format"""
        fieldNew = fieldNew.replace(replaceChar,str(fieldOld)+replaceChar2)
        typeOld = type(fieldOld)
        if isinstance(fieldOld,int) or isinstance(fieldOld,float):
            fieldNew = typeOld(eval(fieldNew))
        else:
            fieldNew = str(fieldNew)
        ###
        return fieldNew

    def isSameName(self):
        """used for nastran = format"""
        if '=' in self.card[0]:
            return True
        return False

    def applyOldFields(self,cardCount=0):
        """used for nastran = format"""
        if not self.isSameName():
           return

        self.cardCount = cardCount
        stop = False
        self.cardTextOld = self.card
        if self.nfields==1 and self.oldCard.nfields>1:
            #self.nfields = self.oldCard.nfields
            #self.applyOldFields(self,cardCount=1)
            cardCount = self.oldCard.cardCount+cardCount
            if self.debug:
                print("newCount = ",cardCount)
            self.card = copy.deepcopy(self.oldCard.cardTextOld)
            self.nfields = len(self.card)
            self.oldCard.nfields = len(self.oldCard.card)
            
            if self.debug:
                print("oldCard = ",self.oldCard)
                print("selfCard = ",self.card)
                print("asfdasdfasdfasdfasfasdf")
            #stop = True

        fieldsNew = self.fields()
        fieldsOld = self.oldCard.fields()
        
        maxLength = max(self.nFields(),self.oldCard.nFields())
        minLength = min(self.nFields(),self.oldCard.nFields())
        
        cardBuilt = [fieldsOld[0]]
        
        for i in range(1,minLength):
            fieldOld = fieldsOld[i]
            fieldNew = fieldsNew[i]

            a = "|%s|" %(fieldNew)
            if '*' in fieldNew:
                newChar = '+%s*'%(cardCount+1)
                fieldNew = self.replaceExpression(fieldNew,fieldOld,'*',newChar)
            elif '/' in fieldNew:
                newChar = '-%s*'%(cardCount+1)
                fieldNew = self.replaceExpression(fieldNew,fieldOld,'/',newChar)
            elif '=='==fieldNew:
                #break
                fieldNew = fieldOld
            elif '=' in fieldNew:
                if fieldNew=='=':
                    fieldNew = fieldOld
                elif fieldNew=='==': # handle this in the max length section
                    pass
                else: # replace = with str(expression)
                    fieldNew = self.replaceExpression(fieldNew,fieldOld)
                ###
            elif ''==fieldNew:
                fieldNew = fieldOld
            else:
                b = "|%s|" %(fieldNew)
                c = "|%s|" %(fieldOld)
                print("i=%s fieldStart %-10s fieldNew %-10s fieldOld %-10s" %(i,a,b,c))
                raise Exception('unhandled case...')
            b = "|%s|" %(fieldNew)
            c = "|%s|" %(fieldOld)
            if self.debug:
                print("i=%s fieldStart %-10s fieldNew %-10s fieldOld %-10s" %(i,a,b,c))
            cardBuilt.append(fieldNew)
            i+=1
        ###
        
        if maxLength<len(cardBuilt): # the new card is longer than builtCard
            for i in range(self.nfields,maxLength):
                cardBuilt.append(self.card[i])
        elif len(cardBuilt)<self.oldCard.nfields: # builtCard is shorter than the old card
            for i in range(self.nfields,maxLength):
                cardBuilt.append(self.oldCard.field(i))
        else: # same length
            pass
        
        if self.debug:
            print("cardBuilt = ",cardBuilt)
        self.card = cardBuilt
        self.nfields = len(self.card)
        
        if stop:
            sys.exit("asdfasdf")

    def getOldField(i):
        """used for nastran = format"""
        return self.oldCard.field(i)

def wipeEmptyFields(card):
    """
    Removes an trailing Nones from the card.  Also converts empty strings to None.
    @param self the object pointer
    @param card the fields on the card as a list
    @retval shortCard the card with no trailing blank fields
    """
    #print "cardA = ",card
    cardB = []
    for field in card:
        if isinstance(field,str) and field.strip()=='':
            field = None
        ###
        cardB.append(field)
    ###
    #print "cardB = ",cardB
    i = 0
    iMax = 0
    while i<len(card):
        if cardB[i] is not None:
            iMax = i
        ###
        i+=1
    ###
    #print "i=%s iMax=%s"%(i,iMax)
    #print "cardC = ",cardB[:iMax+1],'\n'

    return cardB[:iMax+1]
