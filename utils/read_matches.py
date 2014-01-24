

class ReadMatches(object):
    """Class encapsulating the logic of storing matches to reads."""
    def __init__(self, read_id):
        self.read_id = read_id
        
        # Private. Read and written as properties.
        self._insert_match = None
        self._backbone_match = None

        # Insert-related.
        # Strand of the query on which the insert matched.
        self.insert_match_strand = None
        # Range of the query to which insert matched.
        # Ordered as in the BLAT match.
        self.insert_match_range = None
        # End of the insert that was matched. '3p' or '5p'
        self.insert_match_end = None
        # Computed location of the insert in the read.
        self.insert_location_in_query = None

        # Backbone-related.
        # Strand of the query on which the backbone matched.
        self.backbone_match_strand = None
        # Range of the query to which backbone matched.
        # Ordered as in the BLAT match.
        self.backbone_match_range = None
        # Computed location of the backbone in the read.
        self.backbone_location_in_query = None
        self.match_location_in_backbone = None
        
    def Complete(self):
        """Returns true if all data is initialized"""
        return (self._insert_match and
                self._backbone_match and
                self.read_id)

    def Consistent(self):
        """Check that the insert and backbone are reasonably spaced."""
        if not self.Complete():
            return False

        if self.insert_match_strand != self.backbone_match_strand:
            # Matching on opposite strands of query is no good.
            return False

        loc_diff = self.backbone_location_in_query - self.insert_location_in_query
        # If the ends are too far away in the read, remove.
        # TODO: figure out what the correct logic here is.
        return abs(loc_diff) < 6
    
    def get_insert_match(self):
        return self._insert_match

    def set_insert_match(self, match):
        self._insert_match = match
        
        # TODO: unsure if this is the right way to get the strand.
        self.insert_match_strand = match.fragments[0].query_strand
        matched_end = '3p' if match.hit_id.endswith('_3p') else '5p'
        self.insert_match_end = matched_end
        self.insert_match_range = match.query_range
        
        # Correct for off-by-one indexing in BLAT
        corrected_qstart = match.query_start + 1
        if self.insert_match_strand == 1:  # Sense strand.
            loc = corrected_qstart if matched_end == '5p' else match.query_end
        else:
            loc = match.query_end if matched_end == '5p' else corrected_qstart
        self.insert_location_in_query = loc
  
    insert_match = property(get_insert_match, set_insert_match)
    
    def get_backbone_match(self):
        return self._backbone_match
    
    def set_backbone_match(self, match):
        """Set the backbone match.

        Also computes various data from the match including:
            * Strand of match.
            * Location of match in query.
            * Location of match in backbone.
        """
        self._backbone_match = match
        
        # TODO: unsure if this is the right way to get the strand.
        self.backbone_match_strand = match.fragments[0].query_strand
        self.backbone_match_range = match.query_range

        # Correct for off-by-one indexing in BLAT
        insert_end = self.insert_match_end  # Short name
        corrected_qstart = match.query_start + 1
        if self.backbone_match_strand == 1:  # Sense strand
			query_position = match.query_end if insert_end == '5p' else corrected_qstart
        else:
			query_position = corrected_qstart if insert_end == '5p' else match.query_end
        self.backbone_location_in_query = query_position

        corrected_hstart = match.hit_start + 1
        if insert_end == '5p':
            hit_position = match.hit_end - 5
        else:
            # One of the ends needs to be shifted by 5 to account
            # for the transposons copying of 5bp.
            hit_position = corrected_hstart
        self.match_location_in_backbone = hit_position

    backbone_match = property(get_backbone_match, set_backbone_match)

