# Circos Class -------------------------------------------------------------------
Circos <- S7::new_class(
  "Circos",
  properties = list(

    # Writable
    tracks = S7::class_list,
    inter_chrom_spacing = S7::class_numeric,

    # Computed
    n_tracks = S7::new_property(
      class = S7::class_numeric,
      getter = function(self) { length(self@tracks)},
      setter = function(self, value) { "@n_tracks is a read-only property" }
    ),
    tracktypes = S7::new_property(
      class = S7::class_character,
      getter = function(self) {
        vapply(self@tracks, FUN = function(z) { z@type}, FUN.VALUE = character(1))
      },
      setter = function(self, value) { "@tracktypes is a read-only property" }
    ),

    n_reference_tracks = S7::new_property(
      class = S7::class_numeric,
      getter = function(self) { sum(self@tracktypes == "Reference") },
      setter = function(self, value) { "@n_reference_tracks is a read-only property" }
    ),

    reference_names = S7::new_property(
      class = S7::class_list,
      getter = function(self) {
        stdref <- self@tracks[[1]]@data # First track is always reference
        unique(stdref$name)
      },
      setter = function(self, value) { "@reference_names is a read-only property" }
    ),

    # Ref Normalised Coord
    ls_normalised_map = S7::new_property(
      class = S7::class_list,
      getter = function(self) {
        stdref <- self@tracks[[1]]@data
        create_normalisation_map(stdref = stdref, self@inter_chrom_spacing)
        },
      setter = function(self, value) { "@ls_normalised_map is a read-only property" }
    ),

    reference_normalised_data = S7::new_property(
      class = S7::class_data.frame,
      getter = function(self) {
        stdref <- self@tracks[[1]]@data
        chrom <- stdref$name
        start <- stdref$start
        end <- stdref$end

        df = project_intervals(
          chr = chrom,
          start = start,
          end = end,
          map = self@ls_normalised_map,
          return_type = "data.frame"
        )
        df <- dfx::rename(df, namemap = c("name", "chr"))
        return(df)
      },
      setter = function(self, value) { "@n_tracks is a read-only property" }
    ),

    # Functions
    project_points = S7::new_property(
      class = S7::class_function,
      getter = function(self){
        function(chr, position) {
          project_points(chr = chr, pos = position, map = self@ls_normalised_map)
        }
      }
    ),

    project_intervals = S7::new_property(
      class = S7::class_function,
      getter = function(self){
        function(chr, start, end) {
          project_intervals(chr = chr, start = start, end = end, map = self@ls_normalised_map)
        }
      }
    )
  ),
  validator = function(self) {
    if (self@n_reference_tracks != 1) sprintf("Circos class object must have exactly one reference genome, not %i", self@n_reference_tracks)
    if (self@tracktypes[[1]] != "Reference") sprintf("First track must be of type `Reference`, not %s", self@tracktypes[[1]])
    tracks_classes_valid <- all(vapply(self@tracks, function(z) {inherits(z, "ggcircos::Track") }, FUN.VALUE = logical(1)))
    if (!tracks_classes_valid) sprintf("@tracks paramater must be entirely composed of ggcircos::Track objects created by `Track()` function")
  },
  constructor = function(reference, inter_chrom_spacing = 0) {

    # Identify reference
    reference_track <- Track(reference, type = "Reference")

    S7::new_object(
      S7::S7_object(),
      tracks = list(reference_track),
      inter_chrom_spacing = inter_chrom_spacing
    )
})

# Track Class -------------------------------------------------------------------
Track <- S7::new_class("Track",
  properties = list(
    type = S7::new_property(
      class = S7::class_character,
      validator = function(value) {
        if (length(value) != 1) {
          return("@type must be a string (length 1 character vector)")
        }
        if (!value %in% valid_track_types()) {
          return(sprintf("@type must be one of [%s]", toString(valid_track_types())))
        }
      }
    ),
    # A data.frame with columns whose expected cols depend on type
    data = S7::new_property(
      class = S7::class_data.frame
    )
  ),
  constructor = function(data, type) {
    assertions::assert_one_of(type, valid_track_types())
    assertions::assert_dataframe(data)

    if (type == "Reference") {
      assertions::assert_names_include(data, c("name", "start", "end", "refindex"))
    }

    if (type == "SmallMutations") {
      assertions::assert_names_include(data, c("name", "position", "vaf"))
    }

    S7::new_object(
      S7::S7_object(),
      data = data,
      type = type
    )
  }
)


# Enums -------------------------------------------------------------------
valid_track_types <- function() {
  c("Reference", "Copynumber", "SmallMutations", "Links")
}
