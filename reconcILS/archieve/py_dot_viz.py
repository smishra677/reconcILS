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
    