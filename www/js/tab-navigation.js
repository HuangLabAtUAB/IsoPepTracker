/**
 * Tab Navigation JavaScript for IsoPepTracker Shiny Application
 * Originally embedded in ui.R
 * Handles tab switching and active state management
 */

$(document).ready(function() {
  // Initialize active tab
  $('.tab-overview').addClass('active');
  
  // Handle tab clicks
  $('.tab-btn').on('click', function() {
    // Remove active class from all tabs
    $('.tab-btn').removeClass('active');
    
    // Add active class to clicked tab
    $(this).addClass('active');
  });
  
  // Listen for programmatic tab changes
  $(document).on('shiny:inputchanged', function(event) {
    if (event.name === 'canonical_current_tab') {
      // Remove active class from all tabs
      $('.tab-btn').removeClass('active');
      
      // Add active class to the corresponding tab
      if (event.value === 'overview') {
        $('.tab-overview').addClass('active');
      } else if (event.value === 'events') {
        $('.tab-events').addClass('active');
      } else if (event.value === 'isoforms') {
        $('.tab-isoforms').addClass('active');
      }
    } else if (event.name === 'alt_splicing_current_tab') {
      // Remove active class from all tabs
      $('.tab-btn').removeClass('active');
      
      // Add active class to the corresponding tab
      if (event.value === 'rmats') {
        $('.tab-rmats').addClass('active');
      } else if (event.value === 'spladder') {
        $('.tab-spladder').addClass('active');
      }
    }
  });
});