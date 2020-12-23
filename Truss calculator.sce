//Number of truss members and nodes
w_title = "Provide the number of truss members and nodes"
labels = ["Number of members"; "Number of nodes"]
[ok, num_mem, num_nodes] = getvalue(w_title, labels, list("vec", 1, "vec", 1), ["3"; "3"])

//Coordinates of nodes
coord_nodes = zeros(num_nodes,3)
w_title = "Provide the coordinates of node"
labels = ["X coordinate";"Y coordinate"]
for i = 1:num_nodes
    [ok, x_coord, y_coord] = getvalue(w_title + string(i), labels, list("vec", 1, "vec", 1), ["";""])
    coord_nodes(i,:) = [i, x_coord, y_coord]
end

//Member definition
mem_def = zeros(num_mem,3)
w_title = "Provide the definition of member"
labels = ["Start node number"; "End node number"]
for i = 1:num_mem
    [ok, start_m, end_m] = getvalue(w_title + string(i), labels, list("vec", 1, "vec", 1), ["";""])
    mem_def(i,:) = [i, start_m, end_m]
end

//Sketch of the truss
minx = min(coord_nodes(:,2))-1
miny = min(coord_nodes(:,3))-1
maxx = max(coord_nodes(:,2))+1
maxy = max(coord_nodes(:,3))+1

for i=1:num_mem
    xs = coord_nodes(mem_def(i,2),2)
    xe = coord_nodes(mem_def(i,3),2)
    ys = coord_nodes(mem_def(i,2),3)
    ye = coord_nodes(mem_def(i,3),3)
    xv = [xs;xe]
    yv = [ys;ye]
    xsegs(xv,yv)
end

a = gca()
a.data_bounds=[minx,miny,maxx,maxy]

//defining the properties of successive member
Amatrix = zeros(num_mem,1)
Ematrix = zeros(num_mem,1)

w_title = "Provide the properties of member "
labels = ["Cross-section area [cm^2]"; "Young Modulus [GPa]"]
for i = 1:num_mem
    [ok, A, E] = getvalue(w_title + string(i), labels, list("vec", 1, "vec", 1), ["1.396"; "210"])
    Amatrix(i) = A*10^-4
    Ematrix(i) = E*10^9
end

//defining the loads at nodes

loads = zeros(num_nodes, 3)
w_title = "Define the force acting on node "
labels = ["Horizontal component [kN]"; "Vertical component [kN]"]
for i = 1:num_nodes
    [ok, x_comp, y_comp] = getvalue(w_title + string(i), labels, list("vec", 1, "vec", 1), ["0"; "0"])
    loads(i,:) = [i, x_comp*1000, y_comp*1000]
end

//sketching the loads at nodes

max_load = max(abs(loads(:,[2,3])))
for i = 1:num_nodes
    if loads(i,2)<>0
        xv = [coord_nodes(i,2); coord_nodes(i,2)+loads(i,2)/max_load]
        yv = [coord_nodes(i,3); coord_nodes(i,3)]
        xarrows(xv, yv, -1, [3])
    end
    if loads(i,3)<>0
        xv = [coord_nodes(i,2); coord_nodes(i,2)]
        yv = [coord_nodes(i,3); coord_nodes(i,3)+loads(i,3)/max_load]
        xarrows(xv, yv, -1, [3])
    end
end

//defining the reactions

reactions = zeros(num_nodes, 3)
w_title = "Provide the existance of support for node "
labels = ["Horizontal reaction [yes - 1/no - 0]"; "Vertical reaction [yes - 1/no - 0]"]
for i = 1:num_nodes
    [ok, x_comp, y_comp] = getvalue(w_title + string(i), labels, list("vec", 1, "vec", 1), ["0"; "0"])
    reactions(i,:) = [i, x_comp, y_comp]
end

//sketching the reactions at nodes

for i = 1:num_nodes
    if reactions(i,2)<>0
        xv = [coord_nodes(i,2); coord_nodes(i,2)+reactions(i,2)]
        yv = [coord_nodes(i,3); coord_nodes(i,3)]
        xarrows(xv, yv, -1, [2])
    end
    if reactions(i,3)<>0
        xv = [coord_nodes(i,2); coord_nodes(i,2)]
        yv = [coord_nodes(i,3); coord_nodes(i,3)+reactions(i,3)]
        xarrows(xv, yv, -1, [2])
    end
end

//Determining the value of reaction forces

reaction_equations = zeros(3, num_nodes*2)

// Columns go as follows hor_r_node1, hor_r_node2, ..., hor_r_noden, ver_r_node1, ...
for i=1:num_nodes
    reaction_equations(1,i) = reactions(i,2)
    reaction_equations(2,num_nodes+i) = reactions(i,3)
    reaction_equations(3,i) = coord_nodes(i,3)*reactions(i,2)
    reaction_equations(3,num_nodes+i) = -coord_nodes(i,2)*reactions(i,3)
end


sum_of_horizontal_loads = 0
for i=1:num_nodes
    sum_of_horizontal_loads = sum_of_horizontal_loads + loads(i,2)
end

sum_of_vertical_loads = 0
for i=1:num_nodes
    sum_of_vertical_loads = sum_of_vertical_loads + loads(i,3)
end

sum_of_moments = 0
for i=1:num_nodes
    sum_of_moments = sum_of_moments + coord_nodes(i,3)*loads(i,2)-coord_nodes(i,2)*loads(i,3)
end

equilibrium_conditions = [sum_of_horizontal_loads; sum_of_vertical_loads; sum_of_moments]

reaction_values = linsolve(reaction_equations, equilibrium_conditions)

for i=1:num_nodes
    reactions(i,2) = reaction_values(i)
    reactions(i,3) = reaction_values(i+num_nodes)
end

//clearing negligible reactions that should be zero
h = 0.0001
abs(reactions(2,2))<h
for i=1:num_nodes
    for j=1:3
        if abs(reactions(i,j))<0.0001 then
            reactions(i,j) = 0
        end
    end
end

//plotting the new loads that are applied on the truss. We treat reactions as another loading so we add them together

for i=1:num_nodes
    loads(i,2:3)=loads(i,2:3)+reactions(i,2:3)
end

//plotting new truss with combined loading

clf() //we have to clear the current figure because previous reaction arrows are in the way 

minx = min(coord_nodes(:,2))-1
miny = min(coord_nodes(:,3))-1
maxx = max(coord_nodes(:,2))+1
maxy = max(coord_nodes(:,3))+1

for i=1:num_mem
    xs = coord_nodes(mem_def(i,2),2)
    xe = coord_nodes(mem_def(i,3),2)
    ys = coord_nodes(mem_def(i,2),3)
    ye = coord_nodes(mem_def(i,3),3)
    xv = [xs;xe]
    yv = [ys;ye]
    xsegs(xv,yv)
end

a=gca()
a.data_bounds=[minx,miny;maxx,maxy]

max_load = max(abs(loads(:,[2,3])))
for i = 1:num_nodes
    if loads(i,2)<>0
        xv = [coord_nodes(i,2); coord_nodes(i,2)+loads(i,2)/max_load]
        yv = [coord_nodes(i,3); coord_nodes(i,3)]
        xarrows(xv, yv, -1, [3])
    end
    if loads(i,3)<>0
        xv = [coord_nodes(i,2); coord_nodes(i,2)]
        yv = [coord_nodes(i,3); coord_nodes(i,3)+loads(i,3)/max_load]
        xarrows(xv, yv, -1, [3])
    end
end

//calculating internal forces in members

//determining the internal forces 
internal_forces_at_nodes = zeros(num_nodes,1+num_mem*2)
// columns go as follows: N1x, N2x, N3x, ... Nnx, N1y, N2y, ..., Nny

for i = 1:num_nodes
    internal_forces_at_nodes(i,1) = i
end

//defining internal forces acting on each node

 
function L=len(mem)
    nodes = mem_def(mem,2:3)
    node1 = coord_nodes(nodes(1),2:3)
    node2 = coord_nodes(nodes(2),2:3)
    L = sqrt((node1(1)-node2(1))^2+(node1(2)-node2(2))^2)
endfunction


for i=1:num_mem
    node1 = mem_def(i,2)
    node2 = mem_def(i,3)
    coord_node1 = coord_nodes(node1, 2:3)
    coord_node2 = coord_nodes(node2, 2:3)
    internal_forces_at_nodes(node1, 1+i) = (coord_node2(1)-coord_node1(1))/len(i)
    internal_forces_at_nodes(node1, 1+num_mem+i) = (coord_node2(2)-coord_node1(2))/len(i)
    internal_forces_at_nodes(node2, 1+i) = (coord_node1(1)-coord_node2(1))/len(i)
    internal_forces_at_nodes(node2, 1+num_mem+i) = (coord_node1(2)-coord_node2(2))/len(i)
end

//defining the equilibrium equation for each node
equilibrium_equations_each_node = zeros(2*num_nodes,num_mem)
for i = 1:num_nodes
    for j = 1:num_mem
        equilibrium_equations_each_node(i,j)=internal_forces_at_nodes(i,1+j)
    end
    for j = 1:num_mem
        equilibrium_equations_each_node(num_nodes+i,j)=internal_forces_at_nodes(i,1+j+num_mem)
    end
end

equilibrium_condition_at_node = zeros(2*num_nodes,1)
for i = 1:num_nodes
    equilibrium_condition_at_node(i) = loads(i,2)
    equilibrium_condition_at_node(num_nodes+i) = loads(i,3)
end

internal_forces_in_members = linsolve(equilibrium_equations_each_node, equilibrium_condition_at_node)

//clearing negligible internal forces that should be zero
h = 0.0001
for i=1:length(internal_forces_in_members)
    if abs(internal_forces_in_members(i))<h then
        internal_forces_in_members(i) = 0
    end
end
//calculating internal stresses
stress_vector = zeros(num_mem,1)
for i=1:num_mem
    stress_vector(i) = internal_forces_in_members(i)/Amatrix(i)
end

//converting units [N]->[kN], [Pa]->[MPa]
internal_forces_in_members=internal_forces_in_members*0.001
stress_vector=stress_vector*10^-6

//Creating a table for results
results = []
for i= 1:num_mem
    results(i+1)=string(i)
end
results(1,1)="member no."
results(1,2)="Internal force [kN]"
results(1,3)="Internal Stress [MPa]"
results(1,4)="Type of stress"

for i=1:num_mem
    results(i+1,2)=string(abs(internal_forces_in_members(i)))
end
for i=1:num_mem
    if stress_vector(i)>0 then
        results(i+1,4)="Tension"
    elseif stress_vector(i)<0 then
        results(i+1,4)="Compression"
    else
        results(i+1,4)="-"
    end
end

for i=1:num_mem
    results(i+1,3)=string(abs(stress_vector(i)))
end

csvWrite(results,'truss.txt')
