import networkx as nx
from typing import List, Type
from EukMetaSanity.tasks.base.task_class import TaskList

"""
Class takes list of tasks and creates dependency graph of data
Topological sort outputs order in which tasks can be completed

"""


class DependencyGraphGenerationError(BaseException):
    pass


class DependencyGraph:
    def __init__(self, tasks: List[TaskList]):
        self.idx = {task.name: task for task in tasks}
        self.graph = nx.DiGraph()
        self.build_dependency_graph(tasks)
        if not nx.is_directed_acyclic_graph(self.graph):
            raise DependencyGraphGenerationError("There was an error generating dependency graph information for this "
                                                 "pipeline")

    def build_dependency_graph(self, tasks: List[TaskList]):
        self.graph.add_node("root")
        # Add dependencies stored in TaskList object's .requires member
        for task_list in tasks:
            self.graph.add_edge("root", task_list.name)
            for dependency in task_list.requires:
                self.graph.add_edge(dependency, task_list.name)

    def sorted_tasks(self) -> List[Type[TaskList]]:
        sorted_nodes = list(nx.topological_sort(self.graph))
        sorted_nodes.remove("root")
        return [self.idx[name] for name in sorted_nodes]
