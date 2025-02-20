function vertices = findPolyhedronVertices(A, b)
    % Number of variables
    n = size(A, 2);
    % All combinations of constraints
    combs = nchoosek(1:size(A, 1), n);
    % Initialize an empty array to store the vertices
    vertices = [];
    % Solve for each combination of constraints
    for i = 1:size(combs, 1)
        % Select the constraints
        A_eq = A(combs(i, :), :);
        b_eq = b(combs(i, :));
        % Solve the linear system A_eq * x = b_eq
        x = A_eq \ b_eq;
        % Check if the solution is feasible
        if all(A * x <= (b+1e-3))
            % Add the vertex to the list
            vertices = [vertices, x];
        end
    end
    % Remove duplicate vertices
    vertices = unique(vertices', 'rows')';
end