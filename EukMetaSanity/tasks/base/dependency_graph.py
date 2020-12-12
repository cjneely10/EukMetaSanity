import networkx as nx
from collections import namedtuple
from typing import List, Type, Tuple
from EukMetaSanity.tasks.base.task_class import TaskList
from EukMetaSanity.tasks.dependencies import dependencies


Node = namedtuple("Node", ("name", "scope"))


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
            raise DependencyGraphGenerationError(
                "There was an error generating dependency graph information for this pipeline.\n"
                "If you are a user seeing this, this is likely a fatal error - contact the authors of the pipeline\n"
                "If you are a developer, double-check dependency information in your pipeline"
            )

    def _build_dependency_graph(self, tasks: List[TaskList]):
        self.graph.add_node(Node(name="root", scope=""))
        # Add dependencies stored in TaskList object's .requires member
        for task_list in tasks:
            task_node = Node(name=task_list.name, scope=task_list.name)
            self.graph.add_edge(Node(name="root", scope=""), task_node)
            for dependency in task_list.requires:
                self._add_dependency_stack(Node(name=dependency, scope=task_node.name), task_node.name, task_node.name)

    def _add_dependency_stack(self, dependency: Node, main_task_name: str, scope: str):
        self.graph.add_edge(dependency, Node(name=main_task_name, scope=scope))
        dep_task: TaskList = self.idx[dependency.name]
        for dep in dep_task.requires:
            self.graph.add_edge(Node(name=dep, scope=scope), Node(name=dep_task.name, scope=scope))
            self._add_dependency_stack(Node(name=dep, scope=scope), dep_task.name, scope)

    @property
    def sorted_tasks(self) -> List[Tuple[Type[TaskList], str]]:
        """ Run topological sort of all tasks in pipeline and output in order that allows
        for completion of dependencies in required order

        :return: List of TaskList child classes to run
        """
        sorted_nodes: List[Node] = list(nx.topological_sort(self.graph))
        sorted_nodes.remove(Node(name="root", scope=""))
        return [(self.idx[node.name], node.scope) for node in sorted_nodes]
