# Monkey patch for Psych to load RBTree YAML files.
# First function is for Infinity/NaN appearing as numeric values.
# Second function is for RBTree.
require "psych"
require "psych/visitors/visitor"

module Psych
  module Visitors
    class ToRuby < Psych::Visitors::Visitor


      def revive_hash hash, o
        @st[o.anchor] = hash if o.anchor

        o.children.each_slice(2) { |k,v|
          key = accept(k)
          val = accept(v)

          if key == '<<'
            case v
            when Nodes::Alias
              begin
                hash.merge! val
              rescue TypeError
                hash[key] = val
              end
            when Nodes::Sequence
              begin
                h = {}
                val.reverse_each do |value|
                  h.merge! value
                end
                hash.merge! h
              rescue TypeError
                hash[key] = val
              end
            else
              hash[key] = val
            end
          elsif key == "Infinity"
            hash[Float::INFINITY] = val
          elsif key == "NaN"
            hash[Float::NAN] = val
          else
            begin
              hash[key] = val
            rescue => e
              warn("Unrecognized key or value problem")
              binding.pry
            end
          end

        }
        hash
      end

      def visit_Psych_Nodes_Mapping o
        if Psych.load_tags[o.tag]
          return revive(resolve_class(Psych.load_tags[o.tag]), o)
        end
        return revive_hash({}, o) unless o.tag

        case o.tag
        when /^!ruby\/struct:?(.*)?$/
          klass = resolve_class($1) if $1

          if klass
            s = register(o, klass.allocate)

            members = {}
            struct_members = s.members.map { |x| class_loader.symbolize x }
            o.children.each_slice(2) do |k,v|
              member = accept(k)
              value  = accept(v)
              if struct_members.include?(class_loader.symbolize(member))
                s.send("#{member}=", value)
              else
                members[member.to_s.sub(/^@/, '')] = value
              end
            end
            init_with(s, members, o)
          else
            klass = class_loader.struct
            members = o.children.map { |c| accept c }
            h = Hash[*members]
            klass.new(*h.map { |k,v|
              class_loader.symbolize k
            }).new(*h.map { |k,v| v })
          end

        when /^!ruby\/object:?(.*)?$/
          name = $1 || 'Object'

          if name == 'Complex'
            class_loader.complex
            h = Hash[*o.children.map { |c| accept c }]
            register o, Complex(h['real'], h['image'])
          elsif name == 'Rational'
            class_loader.rational
            h = Hash[*o.children.map { |c| accept c }]
            register o, Rational(h['numerator'], h['denominator'])
          elsif name == 'RBTree'
            revive_hash resolve_class($1).new, o
          else
            obj = revive((resolve_class(name) || class_loader.object), o)
            obj
          end

        when /^!(?:str|ruby\/string)(?::(.*))?/, 'tag:yaml.org,2002:str'
          klass   = resolve_class($1)
          members = {}
          string  = nil

          o.children.each_slice(2) do |k,v|
            key   = accept k
            value = accept v

            if key == 'str'
              if klass
                string = klass.allocate.replace value
              else
                string = value
              end
              register(o, string)
            else
              members[key] = value
            end
          end
          init_with(string, members.map { |k,v| [k.to_s.sub(/^@/, ''),v] }, o)
        when /^!ruby\/array:(.*)$/
          klass = resolve_class($1)
          list  = register(o, klass.allocate)

          members = Hash[o.children.map { |c| accept c }.each_slice(2).to_a]
          list.replace members['internal']

          members['ivars'].each do |ivar, v|
            list.instance_variable_set ivar, v
          end
          list

        when '!ruby/range'
          klass = class_loader.range
          h = Hash[*o.children.map { |c| accept c }]
          register o, klass.new(h['begin'], h['end'], h['excl'])

        when /^!ruby\/exception:?(.*)?$/
          h = Hash[*o.children.map { |c| accept c }]

          e = build_exception((resolve_class($1) || class_loader.exception),
                              h.delete('message'))
          init_with(e, h, o)

        when '!set', 'tag:yaml.org,2002:set'
          set = class_loader.psych_set.new
          @st[o.anchor] = set if o.anchor
          o.children.each_slice(2) do |k,v|
            set[accept(k)] = accept(v)
          end
          set

        when /^!map:(.*)$/, /^!ruby\/hash:(.*)$/
          revive_hash resolve_class($1).new, o

        when '!omap', 'tag:yaml.org,2002:omap'
          map = register(o, class_loader.psych_omap.new)
          o.children.each_slice(2) do |l,r|
            map[accept(l)] = accept r
          end
          map

        else
          revive_hash({}, o)
        end
      end
    end
  end
end