# Function to capitalize a string in CMake
include_guard(DIRECTORY)
function(capitalize_string string output_variable)
  string(TOUPPER "${string}" _upper_string)
  string(TOLOWER "${string}" _lower_string)
  string(SUBSTRING "${_upper_string}" 0 1 _start)
  string(SUBSTRING "${_lower_string}" 1 -1 _end)
  set(${output_variable}
      "${_start}${_end}"
      PARENT_SCOPE
      )
endfunction()
