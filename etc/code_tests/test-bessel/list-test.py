parameters = {
  "external_includes": ""
}
parameters["external_includes"] = "test:test2"

includes = set("#include <%s>" % inc for inc in parameters.get("external_includes", ()))

print(includes)

print(set("%s" % s for s in ["test"]))