function group = flattenGroupsAndMoveSensorToTop(propName,groups)
  % Flatten property groups into one group, and move the named
  % property first.
  propertyList = [groups.PropertyList];
  propNameIdx = strcmp(propName,propertyList);
  group = matlab.mixin.util.PropertyGroup([propName,propertyList(~propNameIdx)]);
end
