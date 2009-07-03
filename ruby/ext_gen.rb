# Module that will generate functions to wrap a c++ class.
#
# Vector3f = ExtGen::CClass.new("Vector3f") do
#   define(:method, "double", "Len", [])
#   define(:func, "double", "Lerp", [["Vector3f", "other"], ["double", "t"]])
#   define(:binary_op, "Vector3f", "plus", [["Vector3f", "other"]], "+")
#   define(:assign, "Vector3f", "initialize", [["Vector3f", "other"]])
#   define(:self_method, "void", "Set", [["Vector3f", "other"]])
#   define(:get_attrib, "double", "x", [])
# end
#
# puts Vector3f.gen


module ExtGen

  class CMethod
    def initialize(klass_name, return_type, c_name, rb_name, args)
      @klass_name = klass_name
      @return_type = return_type
      @c_name = c_name
      @rb_name = rb_name
      @args = [[@klass_name, 'self']] + args
    end

    def gen_proto
      "static VALUE #{@klass_name}_#{@c_name}(#{@args.map {|type, name| "VALUE #{name}"}.join(", ")})\n"
    end

    def get_v type, name
      case type
      when "int"
        "\tint #{name}_v = NUM2INT(#{name});"
      when "double"
        "\tdouble #{name}_v = NUM2DBL(#{name});"
      else
        "\tif (TYPE(#{name}) != T_DATA || RDATA(#{name})->dfree != (RUBY_DATA_FUNC)#{type}_free) rb_raise(rb_eRuntimeError, \"#{name} is not a #{type}\");\n" +
        "\t#{type}& #{name}_v = *(#{type}*)DATA_PTR(#{name});"
      end
    end

    def make_v type, name
      case type
      when "double"
        "\tVALUE #{name} = rb_float_new(0);\n" +
        "\t#{type}& #{name}_v = RFLOAT(#{name})->value;"
      else
        "\tVALUE #{name} = #{type}_alloc(c#{type});\n" +
        "\t#{type}& #{name}_v = *(#{type}*)DATA_PTR(#{name});"
      end
    end

    def gen_def
      "\trb_define_method(c#{@klass_name}, \"#{@rb_name}\", RUBY_METHOD_FUNC(#{@klass_name}_#{@c_name}), #{@args.size - 1});\n"
    end

    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
#{make_v @return_type, 'result'}
    result_v = self_v.#{@c_name}(#{@args[1..-1].map{|type, name| "#{name}_v"}.join(", ")});
    return result;
}
CODE
    end

    def gen
      gen_proto + gen_body
    end
  end

  class CMethodFunc < CMethod
    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
#{make_v @return_type, 'result'}
    result_v = #{@c_name}(#{@args.map{|type, name| "#{name}_v"}.join(", ")});
    return result;
}
CODE
    end
  end

  class CMethodBinaryOp < CMethod
    def initialize(klass_name, return_type, c_name, rb_name, args, c_op)
      super(klass_name, return_type, c_name, rb_name, args)
      @c_op = c_op
    end
    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
#{make_v @return_type, 'result'}
    result_v = self_v #{@c_op} #{@args[1][1] + "_v"};
    return result;
}
CODE
    end
  end

  class CMethodUnaryOp < CMethod
    def initialize(klass_name, return_type, c_name, rb_name, args, c_op)
      super(klass_name, return_type, c_name, rb_name, args)
      @c_op = c_op
    end
    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
#{make_v @return_type, 'result'}
    result_v = #{@c_op} self_v;
    return result;
}
CODE
    end
  end

  class CMethodAssign < CMethod
    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
    self_v = #{@args[1][1] + "_v"};
    return self;
}
CODE
    end
  end

  class CMethodSelf < CMethod
    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
    self_v.#{@c_name}(#{@args[1..-1].map{|type, name| "#{name}_v"}.join(", ")});
    return self;
}
CODE
    end
  end

  class CMethodArg1 < CMethod
    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
    self_v.#{@c_name}(#{@args[1..-1].map{|type, name| "#{name}_v"}.join(", ")});
    return #{@args[1][1]};
}
CODE
    end
  end

  class CMethodAttribGet < CMethod
    def initialize *rest
      super
      @attrib_name = @c_name
      @c_name = "get_" + @c_name
    end

    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
#{make_v @return_type, 'result'}
    result_v = self_v.#{@attrib_name};
    return result;
}
CODE
    end
  end

  class CMethodAttribSet < CMethod
    def initialize *rest
      super
      @attrib_name = @c_name
      @c_name = "set_" + @c_name
    end

    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
    self_v.#{@attrib_name} = #{@args[1][1]}_v;
    return #{@args[1][1]};
}
CODE
    end
  end

  class CMethodInitialize < CMethod
    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
#{@args[1..-1].map {|type, name| "\tself_v.#{name} = #{name}_v;"}.join("\n")}
    return self;
}
CODE
    end
  end

  class CMethodCustom < CMethod
    def initialize(klass_name, return_type, c_name, rb_name, args, code)
      super(klass_name, return_type, c_name, rb_name, args)
      @code = code
    end
    def gen_body
<<CODE
{
#{@args.map {|type, name| get_v type, name}.join("\n")}
#{@code}
}
CODE
    end
  end

  class CMethodClass < CMethod
    def gen_def
      "\trb_define_module_function(c#{@klass_name}, \"#{@rb_name}\", RUBY_METHOD_FUNC(#{@klass_name}_#{@c_name}), #{@args.size - 1});\n"
    end

    def gen_body
<<CODE
{
#{@args[1..-1].map {|type, name| get_v type, name}.join("\n")}
#{make_v @return_type, 'result'}
    result_v = #{@klass_name}::#{@c_name}(#{@args[1..-1].map{|type, name| "#{name}_v"}.join(", ")});
    return result;
}
CODE
    end
  end

  MethodFactory = {:method => CMethod, 
                   :self_method => CMethodSelf,
                   :arg1_method => CMethodArg1,
                   :func => CMethodFunc, 
                   :binary_op => CMethodBinaryOp, 
                   :unary_op => CMethodUnaryOp,
                   :assign => CMethodAssign,
                   :get_attrib => CMethodAttribGet,
                   :set_attrib => CMethodAttribSet,
                   :initialize => CMethodInitialize,
                   :custom => CMethodCustom,
                   :class_method => CMethodClass}

  class CClass
    def initialize(klass_name, &block)
      @klass_name = klass_name
      @methods = []
      instance_eval &block
    end

    def define(kind, return_type, c_name, rb_name, args, *rest)
      @methods.push(MethodFactory[kind].new(@klass_name, return_type, c_name, rb_name, args, *rest))
    end

    def gen_functions
memory_functions = <<CODE
static void #{@klass_name}_free(void* p)
{
	delete reinterpret_cast<#{@klass_name}*>(p);
}

// alloc
static VALUE #{@klass_name}_alloc(VALUE klass)
{
	#{@klass_name}* vec = new #{@klass_name}();
	return Data_Wrap_Struct(klass, 0, #{@klass_name}_free, vec);
}
CODE
      memory_functions + "static VALUE c#{@klass_name};\n" + @methods.map {|m| m.gen}.join()
    end

    def gen_defs
      "\tc#{@klass_name} = rb_define_class(\"#{@klass_name}\", rb_cObject);\n" +
      "\trb_define_alloc_func(c#{@klass_name}, #{@klass_name}_alloc);\n" +
      @methods.map {|m| m.gen_def}.join()
    end

  end

end
