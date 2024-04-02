import pickle
import sys
sys.path.append("../")
import utils.Tree as Tree

from pympler import asizeof



# Create an instance of MyClass
original_instance = Tree.Tree()
original_instance.id=2


serialized_instance = pickle.dumps(original_instance)
deep_copied_instance = pickle.loads(serialized_instance)

print(deep_copied_instance,original_instance,deep_copied_instance.id)

print("Size using pympler.asizeof.asizeof():", asizeof.asizeof(deep_copied_instance), "bytes")