function target_vertex = bfs_distance(mesh, start_vertex, distance)
% mesh is a struct that represents a triangulated mesh
% start_vertex is the index of the starting vertex
% distance is the desired distance from the start vertex to the target vertex

visited = false(size(mesh.vertices, 1), 1);
visited(start_vertex) = true;

queue = [start_vertex, 0];

while ~isempty(queue)
    vertex = queue(1, 1);
    dist = queue(1, 2);
    queue(1, :) = [];

    if dist == distance
        target_vertex = vertex;
        return;
    end

    % Find the neighboring vertices
    neighbors = unique(mesh.faces(any(mesh.faces == vertex, 2), :));
    neighbors(neighbors == vertex) = [];

    for i = 1:length(neighbors)
        neighbor = neighbors(i);
        if ~visited(neighbor)
            visited(neighbor) = true;
            queue = [queue; neighbor, dist+1];
        end
    end
end

% If no vertex was found at the desired distance, return -1
target_vertex = -1;
end

