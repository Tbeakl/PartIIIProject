from typing import Type

class Node:

    def __init__(self, name):
        self.name = name
        self.neighbours = []

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.name}, [{', '.join([n.name for n in self.neighbours])}])"
    
    def is_valid_neighbour(self, neighbour) -> bool:
        raise NotImplemented()
    
    def add_neighbour(self, neighbour) -> None:
        assert self.is_valid_neighbour(neighbour)
        self.neighbours.append(neighbour)

class Factor(Node):
    
    def __init__(self, name):
        super().__init__(name)
        self.data = None

    def is_valid_neighbour(self, neighbour) -> bool:
        return isinstance(neighbour, Variable)
    

class Variable(Node):
    def is_valid_neighbour(self, neighbour) -> bool:
        return isinstance(neighbour, Factor)

