from typing import ClassVar, Dict, List

class BDFCard:
    x: ClassVar[int] = 0  # Class variable only

class Coord:
    type: ClassVar[str] = ''

class Node:
    type: ClassVar[str] = ''
    cp: ClassVar[int] = 0

class Element:
    type: ClassVar[str] = ''
    nodes: ClassVar[List[int]] = [1]

class Property:
    type: ClassVar[str] = ''
    pid: ClassVar[int] = 1

class Material:
    type: ClassVar[str] = ''
    mid: ClassVar[int] = 1

class BDF:
    bdf_filename: ClassVar[str] = ''
    nodes : ClassVar[Dict[int, Node]] = {}
    elements : ClassVar[Dict[int, Element]] = {}
    properties : ClassVar[Dict[int, Property]] = {}
    materials : ClassVar[Dict[int, Material]] = {}

