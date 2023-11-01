import igraph as ig

import matplotlib.pyplot as plt


g = ig.Graph()

g = ig.Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)])

g.vs["name"] = ["Alice", "Bob", "Claire", "Dennis", "Esther", "Frank", "George"]
g.vs["age"] = [25, 31, 18, 47, 22, 23, 50]
g.vs["gender"] = ["f", "m", "f", "m", "f", "m", "m"]
g.es["is_formal"] = [False, False, True, True, True, False, True, False, False]
g.vs["label"] = g.vs["name"]

fig, ax = plt.subplots()
layout = g.layout("kk")
ig.plot(g, layout=layout, target=ax)
plt.show()

def visualize_binary_tree(self, graph):
        
        if self.leftChild:
            if len(str(self.leftChild.taxa))==0:
                label='None'
            else:
                label=str(self.leftChild.taxa)
            Node_left=pydot.Node(str(self.leftChild.id), label=label)
            node2 = graph.add_node(Node_left)
            if self.leftChild.evolve=='Loss':
                graph.add_edge(pydot.Edge(str(self.id),str(self.leftChild.id), color='blue'))
            else:
                graph.add_edge(pydot.Edge(str(self.id), str(self.leftChild.id), color='green'))
            self.leftChild.visualize_binary_tree(graph)
        if self.rightChild:
            if len(str(self.rightChild.taxa))==0:
                label='None'
            else:
                label=str(self.rightChild.taxa)
            Node_right=pydot.Node(str(self.rightChild.id), label=str(self.rightChild.taxa))
            node3 = graph.add_node(Node_right)
            if self.rightChild.evolve=='Loss':
                graph.add_edge(pydot.Edge(str(self.id),str(self.rightChild.id), color='blue'))
            else:
                graph.add_edge(pydot.Edge(str(self.id), str(self.rightChild.id), color='green'))
            self.rightChild.visualize_binary_tree(graph)

        else:
            return graph
        
        

        
        
        return graph


    def viz(self):
        graph = pydot.Dot()
        Node_root=pydot.Node(str(self.id), label=('None'))
        node1 = graph.add_node(Node_root)
        #graph.write_png("binary_tree.png")
        graph= self.visualize_binary_tree(graph)
        print(graph)
        graph.write_png("binary_tree.png")
    


def find_all_edges(tree):
  edges = []
  for node in tree:
    for child in node.children:
      edges.append((node, child))
  return edges