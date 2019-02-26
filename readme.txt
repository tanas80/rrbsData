After adding new @export(ed) functions call
call devtools::document()
or roxygen2::roxygenize() # updates .Rd file in /man
to update file NAMESPACE with new exported entries.
