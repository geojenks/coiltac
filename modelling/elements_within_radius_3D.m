function result = elements_within_radius_3D(matrix, center, radius)
    % Get the size of the matrix
    [rows, cols, aisles] = size(matrix);
    
    % Initialize an empty array to store the result
    result = [];
    
    % Loop through each element of the matrix
    for i = 1:rows
        for j = 1:cols
            for k = 1:aisles
                % Calculate the distance between the current element and the center
                distance = sqrt((i - center(1))^2 + (j - center(2))^2 + (k - center(3))^2);
                
                % Check if the distance is less than or equal to the radius
                if distance >= radius(1) && distance <= radius(2) 
                    % Add the element to the result array
                    result = [result, matrix(i, j, k)];
                end
            end
        end
    end
end