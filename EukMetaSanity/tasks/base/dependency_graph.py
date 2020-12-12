import networkx as nx
from typing import List, Type
from EukMetaSanity.tasks.base.task_class import TaskList
from EukMetaSanity.tasks.dependencies import dependencies


class DependencyGraphGenerationError(BaseException):
    """ Exception class if unable to form directed acyclic (dependency) graph from
    provided pipeline configuration

    """
    pass


class DependencyGraph:
    """ Class takes list of tasks and creates dependency graph of data
    Topological sort outputs order in which tasks can be completed

    All tasks are assumed to require (at least) a root task to complete

    The root task is simply populating input files to run in a pipeline

    """
    def __init__(self, tasks: List[TaskList]):
        self.idx = {task.name: task for task in tasks}
        self.idx.update(dependencies)
        self.graph = nx.DiGraph()
        self._build_dependency_graph(tasks)
        if not nx.is_directed_acyclic_graph(self.graph):
            raise DependencyGraphGenerationError("There was an error generating dependency graph information for this "
                                                 "pipeline")

    def _build_dependency_graph(self, tasks: List[TaskList]):
        self.graph.add_node("root")
        # Add dependencies stored in TaskList object's .requires member
        for task_list in tasks:
            self.graph.add_edge("root", task_list.name)
            for dependency in task_list.requires:
                self.graph.add_edge(dependency, task_list.name)

    def sorted_tasks(self) -> List[Type[TaskList]]:
        """ Run topological sort of all tasks in pipeline and output in order that allows
        for completion of dependencies in required order

        :return: List of TaskList child classes to run
        """
        sorted_nodes = list(nx.topological_sort(self.graph))
        sorted_nodes.remove("root")
        return [self.idx[name] for name in sorted_nodes]
